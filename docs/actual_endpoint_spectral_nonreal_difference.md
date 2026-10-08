# Actual endpoint pairs with nonreal spectral differences

Endpoint closeness, primitive Gaussian factors, and comparable logarithmic
factor sizes do not force the displayed factors' slope polynomials to
share leading moments or have real-rooted differences. The following
exact family has three growing factors, pairwise coprime rational norms,
and distinct positive limiting slopes. Its two lattice points have
endpoint arc constant tending to zero, while the difference of their
slope polynomials has two nonreal roots tending to `+i,-i`.

This is a counterfixture to that specific transfer. There are only two
rows: all three factors represent the same unoriented row cut, and
canonical cut grouping combines them into one factor. Equal logarithmic
weights of these three factors are not a near-uniform eight-row cut
profile. Thus this does not obstruct a transfer using a different
grouping, prove anything against the eight-row extraction in
[uniform_profile_extraction.md](uniform_profile_extraction.md), or rule
out other spectral methods.

## An exact primitive factor family

For `k>=1`, define positive integers `x=x_k,y=y_k` by

```
y+x sqrt(75) = (1351+156 sqrt(75))^k.
```

The norm of the displayed unit is one, so

```
y^2-75x^2=1,
(x',y')=(1351x+156y,11700x+1351y).
```

Starting at `(x_0,y_0)=(0,1)`, these recurrences give

```
x even, y odd, gcd(x,y)=1,
y=1 (mod 5), y=(-1)^k (mod 13).
```

In particular `x>=156`, and `sqrt(75)<y/x<9`. Put

```
Gamma_- = x+i(y-8x),
Gamma_0 = 5x+iy,
Gamma_+ = x+i(y+8x).
```

All real and imaginary parts are positive. Each factor has coprime
coordinates: this is immediate for `Gamma_-,Gamma_+`, and for `Gamma_0`
also uses `5` not dividing `y`. Each norm is odd. Consequently each
factor is coprime in `Z[i]` to its conjugate. Direct multiplication gives

```
Gamma_- Gamma_+ = -10x^2-1+2ixy,
G := Gamma_- Gamma_0 Gamma_+ = -x(200x^2+7)-iy.       (1)
```

Write `X=x(200x^2+7)`. The product is also conjugate-primitive.
Indeed, `gcd(x,y)=1` and

```
75(200x^2+7)=200y^2+325.
```

Every prime common to `X,y` would therefore divide `325`, whereas
neither `5` nor `13` divides `y`. Thus `gcd(X,y)=1`; its norm is odd.

The rational norms are

```
N_- = 140x^2-16xy+1,
N_0 = 100x^2+1,
N_+ = 140x^2+16xy+1.                                (2)
```

Here is a direct proof that they are pairwise coprime, including
possible conjugate Gaussian prime factors. If a rational prime `p`
divides both `N_0` and `N_±`, then `p!=2,5`. In `F_p`,

```
x^2=-1/100, y^2=1/4, 16xy=+2/5 or -2/5.
```

Squaring the last equality gives `-16/25=4/25`, impossible for a prime
other than `2,5`. If `p` divides `N_-` and `N_+`, their oddness and
difference `32xy` imply `p|xy`. Their residues modulo any prime divisor
of `x` are one, so `p|y`. The Pell relation then makes
`N_±=65x^2 (mod p)`, and hence `p=5` or `13`, again impossible.

All three norms grow at the same logarithmic rate:

```
N_-/x^2 -> 140-80 sqrt(3)>0,
N_0/x^2 -> 100,
N_+/x^2 -> 140+80 sqrt(3),
log N_j / log(N_- N_0 N_+) -> 1/3.                  (3)
```

## Actual points and the endpoint constant

Take the two integer points

```
z_0=-G=X+iy, z_1=-bar(G)=X-iy,
R=|z_0|=sqrt(X^2+y^2).
```

They have no common Gaussian divisor and have the displayed common
Gaussian unit `-1` in the block factorization. Their smaller angular
separation and their arc length divided by `sqrt(R)` are

```
Delta=2 arctan(y/X),
C_k=Delta sqrt(R).
```

As `k` tends to infinity,

```
R/(x^3) -> 200,
x^2 Delta -> sqrt(75)/100,
sqrt(x) C_k -> sqrt(3/2).                           (4)
```

In particular `C_k->0`. An elementary bound, avoiding an asymptotic
estimate, is `C_k <= (27/20)/sqrt(x)`: use `y<=9x`, `X>=200x^3`,
`R<=X+y<=225x^3`, and `arctan(t)<=t` for positive `t`.

The common-unit statement can also be obtained at the split-prime level
if Gaussian prime generators may be chosen, as in the usual allocation
representation. The integer `N_0=(10x)^2+1` lies strictly between
consecutive squares. Its prime factorization therefore has an odd
exponent; pairwise norm coprimality preserves that odd exponent in
`Norm(G)`. Since `G` is conjugate-primitive with odd norm, this is an odd
exponent of a split Gaussian prime. Adjusting the chosen generator of
that prime by a unit absorbs any total unit in the factorization of
`G`. Then `G` and `bar(G)` are conjugate products with unit one, and
both displayed source points have unit `-1`. No assertion about a
preassigned canonical choice of prime generators is needed.

## Exact spectral identity and failure of moment transfer

Let `q=y/x` and take the literal imaginary-to-real slopes of the three
displayed factors:

```
a_- = q-8, a_0=q/5, a_+=q+8.
```

They are strictly positive and increasing, since `8<q<9<10`, and their
limits are the three distinct numbers

```
5 sqrt(3)-8, sqrt(3), 5 sqrt(3)+8.
```

Define the two monic real-rooted polynomials

```
P_+(t)=product_j(t-a_j),
P_-(t)=product_j(t+a_j)=-P_+(-t).
```

They have the common spectral product

```
P_+(t)P_+(-t)=P_-(t)P_-(-t)=-product_j(t^2-a_j^2).
```

Nevertheless, with `e_r` the elementary symmetric functions of the
three positive slopes,

```
e_1=11q/5,
e_3=q(11+1/x^2)/5,
P_+(t)-P_-(t)=-2(e_1 t^2+e_3),
e_3/e_1=1+1/(11x^2).                               (5)
```

The difference has exactly the two nonreal roots

```
t=+i sqrt(1+1/(11x^2)), -i sqrt(1+1/(11x^2)).        (6)
```

The first odd moments of the two root lists are `+e_1,-e_1`, whose
difference tends to `22 sqrt(3)`, rather than zero. At the complex
evaluation point encoding the actual Gaussian product, however,

```
P_+(i)=-y/(5x^3)+iX/(5x^3),
P_-(i)=+y/(5x^3)+iX/(5x^3),
P_+(i)-P_-(i)=-2y/(5x^3) -> 0.                     (7)
```

Thus endpoint proximity here corresponds to nonreal difference roots
approaching the evaluation points `+i,-i`. The common spectral product
and uniformly separated limiting slopes do not turn this small
evaluation into shared leading coefficients or real-rooted differences.
Those additional hypotheses in
[critical_moment_spectral_orientation.md](critical_moment_spectral_orientation.md)
remain substantive; the present family does not satisfy them.

The exact checker
[check_actual_endpoint_spectral_nonreal_difference.py](check_actual_endpoint_spectral_nonreal_difference.py)
tests 32 Pell iterates using integer and rational arithmetic, including
the literal Gaussian products, rational norm coprimality, polynomial
products, complex evaluations, and the exact nonreal-root formula.

## A direct Herglotz kernel does not transfer endpoint closeness

There is an elementary positive-kernel description at the actual
evaluation point, without imposing moment identities. For any real
slopes `a_c` and sign rows `S_i`, put

```
P_i(t)=product_c(t-S_ic a_c),
H_i(t)=-P_i'(t)/P_i(t).
```

Each `H_i` maps the upper half-plane into itself, since it is a sum of
terms `-1/(t-r)` with real `r`. Exactly at `i`,

```
H_i(i)=A_i+iB,
A_i=sum_c S_ic a_c/(1+a_c^2),
B=sum_c 1/(1+a_c^2).                                (8)
```

All imaginary parts coincide. For actual block coordinates
`Gamma_c=x_c+i y_c`, `a_c=y_c/x_c`, the source normalization is
`z_i=epsilon i^(-m) (product_c x_c) P_i(i)`, whenever the displayed
block factorization has common unit `epsilon`. Thus the endpoint
condition controls differences of the values `P_i(i)`, not their
logarithmic derivatives in (8).

The preceding Pell family makes that distinction quantitative. Its two
Herglotz functions have values `+A(q)+iB(q),-A(q)+iB(q)`, where

```
D(q)=q^4-126q^2+4225,
A(q)=2q(q^2-63)/D(q)+5q/(q^2+25),
B(q)=2(q^2+65)/D(q)+25/(q^2+25).
```

Consequently

```
A(q) -> 11 sqrt(3)/20,
B(q) -> 19/20,
|H_+(i)-H_-(i)| -> 11 sqrt(3)/10.                    (9)
```

Meanwhile the normalized source separation is
`|z_0/R-z_1/R|=2y/R asymp R^(-2/3)`. In particular even separation
smaller than the endpoint scale `R^(-1/2)` does not force the two
logarithmic derivatives to approach one another. This remains a
statement about the displayed refinement, with the scope explained
at the start of the note.

The natural cross-row Cauchy Gram matrix also has a simple exact form.
Set `f_ic=1/(i-S_ic a_c)` and `K_ij=sum_c f_ic conjugate(f_jc)`.
Then `K` is positive semidefinite and

```
K_ij=sum_c [1+a_c^2 S_ic S_jc+i a_c(S_jc-S_ic)]/(1+a_c^2)^2.
```

For any real matrix `L` satisfying `L 1=0`, this reduces to

```
L K L^t = L S diag(w_c) S^t L^t,
w_c=a_c^2/(1+a_c^2)^2=x_c^2 y_c^2/N_c^2,
N_c=x_c^2+y_c^2.                                    (10)
```

Thus this particular centered kernel depends only on factorwise
weighted cut signs. Its positivity holds before imposing endpoint
closeness. For nonzero integer `x_c,y_c`, the elementary bounds are

```
(N_c-1)/N_c^2 <= w_c <= 1/4.
```

In a primitive block representation `product_c N_c=R^2`, a
near-uniform profile with `b` canonical cuts gives only
`w_c >= R^(-2/b+o(1))` from this estimate, while weights may also
remain bounded away from zero. The matrix identity supplies no
endpoint-scale upper bound on these weights or on the derivatives.
A useful continuation would need an additional relation between the
clustered values and this derivative kernel that survives (9), using
many canonical rows and their arithmetic; such a relation is not
established here. This is an audit of one explicit kernel, not a
claim that all positive-kernel methods reduce to (10).
