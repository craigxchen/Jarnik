# Shifted and scaled six-node central ratios

This note audits two simple rational families for six circle nodes. It is
an exact finite or symbolic calculation, not an asymptotic statement about
all cotangent cliques.

## 1. Closed central-weight formula

Anchor the circle at `r_0=1`, and write the five positive cotangent ratios

```text
r_i=(t_i-i)/(t_i+i),       1<=i<=5.
```

The central column is proportional over `Q` to

```text
c_0=1,
c_i=-(t_i^2+1)^2 / product_(j!=i)(t_i-t_j).             (1)
```

Indeed `r_i-r_j=2i(t_i-t_j)/((t_i+i)(t_j+i))`, while
`1-r_i=2i/(t_i+i)`; substituting these into
`r_i^2/product_(j!=i)(r_i-r_j)` gives (1). Clearing the rational
denominators and dividing by their common integer content gives the
primitive ordinary central vector.

## 2. Translates cannot make the ratio small

First let `A` run through positive integers and take `t_i=A+rho_i`,
with fixed distinct rational offsets `rho_i`.
The denominators in (1) are the fixed numbers

```text
D_i=product_(j!=i)(rho_i-rho_j).
```

After multiplying by one fixed integer clearing all polynomial coefficient
denominators (including those of the `rho_i`), the anchor
coordinate stays a fixed nonzero integer, whereas every finite coordinate
is a degree-four polynomial in `A` with nonzero leading coefficient
`-1/D_i`. The common content can divide the fixed anchor coordinate only.
Consequently

```text
max|u_i| / min t_i  grows on the order of A^3.             (2)
```

Thus a fixed-offset translation family cannot have
`max|u_i|/min t_i -> 0`, even before imposing collision avoidance.

## 3. Homotheties have the same obstruction

For positive integer `A` and `t_i=A rho_i`, with fixed distinct positive rational `rho_i`,
clear the fixed denominators in `rho_i` and
`D_i=product_(j!=i)(rho_i-rho_j)`. Formula (1) becomes, up to one fixed
scalar,

```text
c_0=A^4,
c_i=-(A^4 rho_i^4+2A^2 rho_i^2+1)/D_i.                (3)
```

Any common divisor carrying an unbounded factor of `A` would have to
divide the constant `1/D_i` term for every `i`; after the fixed
denominator clearing, its `A`-dependent part is therefore bounded by a
constant depending only on the `rho_i`. Hence the primitive vector still
has size `>> A^4`, while `min t_i` is comparable to `A`, giving the same
`A^3` growth. The observed large ratios in scaled consecutive families
are this explicit lower-order-term obstruction.

The common-content assertion can be stated precisely by a fixed resultant:
the cleared anchor is a nonzero constant times `A^4`, and any one cleared
finite coordinate has a nonzero constant term. These two polynomials are
coprime over `Q[A]`. An integer polynomial Bezout identity bounds their
value gcd by a fixed nonzero integer, and hence bounds the full row gcd.

## 4. Arbitrary rational parameters retain the obstruction

Write `A=s/t` in lowest terms, `t>0`, and consider `A->+infinity`.
After the same fixed coefficient clearing and multiplication by `t^4`,
both families give homogeneous integer quartics in `(s,t)`.
For translation, the anchor is `B t^4` and a finite row is a fixed
nonzero multiple of
`((s+rho_i t)^2+t^2)^2`. For homothety, the anchor is `B s^4` and a
finite row is a fixed nonzero multiple of `(rho_i^2 s^2+t^2)^2`.

In each case the anchor and any one finite row have no common projective
zero over the algebraic closure. Homogeneous Bezout identities give a
fixed nonzero integer `E` such that their common value divisor at every
primitive `(s,t)` divides `E`. Their maximum absolute value on the real
parameter set `max(|s|,|t|)=1` has a positive minimum. Thus primitive row height
is comparable to `max(|s|,|t|)^4`, with constants depending on the fixed
offsets or scale ratios. Along `A->infinity` this gives

```text
max|u_i| / min t_i  asymp s^3 t = A^3 t^4.
```

In particular the ratio is at least a fixed multiple of `A^3`, even
when the denominator of `A` grows. The assertion of exact `A^3` growth
in Sections 2 and 3 is for integer `A`; arbitrary rational parameters
can have a larger ratio.

## 5. Exact bounded scan

The companion checker
[check_positive_six_ratio_search.py](check_positive_six_ratio_search.py)
scans all `142506=binom(30,5)` integer tuples
`1<=t_1<...<t_5<=30`. The unrestricted minimum of
`max|u_i|/t_1` is

```text
289/3 at t=(3,4,8,13,21),
```

but this tuple has exact two-factor product collisions. Every tuple with
ratio below `841/2` has a collision after full primitive Gaussian
normalization. The first collision-free tuple is

```text
t=(2,3,5,12,17),   L=126,
X=(252,378,630,1512,2142),
u=(126,-7,50,-169,841,-841),
N=R^2=273325=5^2*13*29^2,
max|u_i|/min(X_i/L)=841/2.
```

All 15 pair products are distinct even up to Gaussian units. This finite
scan supplies no universal lower bound. The previously used elliptic
fixed-weight family still has pair-product collisions after positive
reanchoring, so it does not provide the requested collision-free
asymptotic family.
