# A direct radial slice with improved logarithmic growth

Let `m` distinct lattice points lie on a circle of radius `R=sqrt(N)`
and on an arc of length at most `C sqrt(R)`, for fixed `C>0`. Choose
any actual point `z_0=x_0+i y_0` on the arc and write

```text
g=gcd(|x_0|,|y_0|),       z_0=g v,
Q=Norm(v)=N/g^2,         rho=sqrt(Q)=R/g,
W=log N,                D=floor(C^2 rho/2).
```

Here `rho` is the length of the primitive integer vector in the anchor's
radial direction. No Gaussian-primitivity, parity, or prime-support
assumption on the full tuple is needed.

For `N>1` and

```text
rho <= log N/(50 C^2),                              (1)
```

we prove

```text
m <= 2 pi(D)+3       if D>=16,
m <= 2D+1 <=31       if D<16,                        (2)
```

where `pi(D)` is the number of rational primes at most `D`. In particular,
throughout (1),

```text
m = O_C(1+rho/log(2+rho)).                           (3)
```

This is a bound for a specified family of arcs. An arbitrary arc need
not contain an anchor satisfying (1), and (3) still grows when `rho`
grows. It does not prove the radius-independent endpoint bound.

## 1. An exact integer model at any actual anchor

For each nonanchor point `z_i`, put

```text
s_i=(N-Re(conj(z_0) z_i))/g=N/g-Re(conj(v) z_i),
t_i=Im(conj(v) z_i),       Z=2N/g=2sqrt(N) rho.
```

Since `g^2|N`, all three quantities are integers. Taking the norm of
`conj(v)z_i` gives the exact equation

```text
t_i^2 = QN-(N/g-s_i)^2 = s_i(Z-s_i).                (4)
```

Also

```text
0<s_i=|z_i-z_0|^2/(2g)<=C^2 R/(2g)=C^2 rho/2.      (5)
```

Every fixed `s_i` gives at most two values of `t_i`, hence at most two
physical points: multiplication by the nonzero Gaussian integer
`conj(v)` is injective. Thus the number of distinct positive deficits
is at least `ceil((m-1)/2)`. If `D=0`, there are no nonanchor points.

## 2. A finite alternative from squareclasses

Assume `D>=1`, and set `r=pi(D)`. Encode each distinct deficit `s` by
the parities of its prime exponents on the `r` primes at most `D`,
followed by one coordinate equal to one. These vectors belong to
`F_2^(r+1)`.

If there are more than `r+1` distinct deficits, a subset of at most
`r+2` vectors has a nonempty linear dependence. The last coordinate
forces its size to be even, say `2k<=r+2`; the other coordinates force
the product of its deficits to be an integer square. Multiplication of
(4) over this subset gives

```text
P(Z)=product_(j=1)^(2k)(Z-s_j)
    =(product_j t_j/sqrt(product_j s_j))^2.         (6)
```

The right side is a rational square and the left side is an integer,
so it is an integer square. The roots are distinct positive integers
at most `D`. If `Z>D`, the
[elementary even-product bound](runge_even_product_square_bound.md)
gives `Z<=(8kD)^(2k+1)`. If `Z<=D`, the following weaker bound is
already automatic. In either case we have the finite alternative

```text
either m<=2r+3,
or     Z<=(8(r+2)D)^(r+3).                         (7)
```

The `+3` counts the actual anchor and the possible two points per
deficit. No independence of the original conductor labels is assumed.

## 3. An elementary prime bound and an explicit radial cutoff

Write `theta(x)=sum_(p<=x)log p`. Every prime in `(n,2n]` divides
`binom(2n,n)`, whose logarithm is at most `2n log 2`. Summing this
inequality over dyadic integer intervals and enclosing `x` in the next
power of two proves

```text
theta(x)<=4(log 2)x.
```

For `x>1`, split the prime count at `sqrt(x)` to get

```text
pi(x)<=sqrt(x)+2theta(x)/log x
     <=(2+8log 2)x/log x <8x/log x,                 (8)
```

using `log x<=2sqrt(x)`. For `D>=16`, (8) implies

```text
r+3<=11D/log D,
log(8(r+2)D)<=4log D.
```

For the second inequality, even the trivial `r<=D` suffices:
`8(D+2)D<=D^4` for `D>=16`. The first follows from (8) and
`3log D<=3D`. Therefore the second branch of (7) implies

```text
log Z<=44D.                                        (9)
```

Under (1), `D<=W/100`, whereas the exact expression for `Z` gives
`log Z=log 2+W/2+log rho>=W/2`. Thus (9) is impossible, since
`44D<=0.44W<W/2`. This proves the first part of (2). For `D<16`,
directly counting the possible integer deficits gives its second part.
Equation (3) follows from (8), with bounded `D` absorbed in the constant.

## 4. Growth consequences and the remaining gap

For every fixed `C`, the condition `rho<=log R/loglog R` eventually
implies (1), and (3) gives

```text
m=O_C(log R/(loglog R)^2).                          (10)
```

More generally, for any fixed `0<sigma<1`, the range
`rho<=(log R)^sigma` gives

```text
m=O_(C,sigma)(1+(log R)^sigma/loglog R).             (11)
```

The earlier [effective radial slice](endpoint_radial_slice_audit.md)
is stronger in its smaller permitted range: it proves a uniform count
using simultaneous Pell equations. Equations (10)--(11) extend improved
growth estimates to larger radial lengths; they do not improve that
earlier uniform count.

Unlike the [nonconformal-map growth theorem](reciprocal_squareclass_runge_growth.md),
this construction uses only the original physical circle. It also
differs from the
[small-content affine obstruction](direct_projection_affine_normalization_height_tradeoff.md):
here a large actual anchor gcd is precisely what supplies the small
deficits. No argument has shown that every large endpoint tuple contains
such an anchor. The general bound and the uniform endpoint goal remain
unchanged.

Root, Astra and Luna independently checked the integer normalization
and the finite squareclass alternative. The
[arithmetic checker](check_anchor_radial_deficit_squareclass_growth.py)
checks the dictionary and square-product passage on literal circles,
including even norms and tuples with a nontrivial Gaussian common factor.
It also retains the five-point Pell control

```text
b^2-2a^2=1,       g=a^2+1,
(g,0), (g-1,b), (g-1,-b), (g-2,2a), (g-2,-2a).
```

All five points have norm `g^2`. At the first anchor the only positive
deficits are `1,2`, each repeated twice. Their augmented parity vectors
are independent, so this family does not supply a square dependence.
The proof removes repeated deficits before applying the polynomial lemma.
