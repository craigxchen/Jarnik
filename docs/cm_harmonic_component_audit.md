# Audit of the inverse harmonic component

This note records an exact algebra check of the inverse parametrization in
`cm_disjoint_mixed_pencil_search.md`.  Write

```text
A=(x-l)^2,  B=(x-m)^2,
P0=x^2,    Q0=x-1,
Pb=c A,    Qb=c(A-r B).
```

The scale `c` does not change `fb=Pb/Qb`.  For the resultant calculation it
is enough to set `c=1`; the projective pencil is unchanged.  In the affine
parameter convention, this gives

```text
Pk=(1-k)P0+k Pb,  Qk=(1-k)Q0+k Qb.
```

Let `Dy` be the discriminant in `x` of `Pk-y Qk`.  Exact expansion gives

```text
D4 = k (a1+a2 k),
D-2 = 12+b1 k+b2 k^2,
```

where

```text
a1=12(l-2)^2+16r(4m-m^2-4)
a2=-12l^2+48l-48+r(48l^2-96lm+64m^2+64-64m)
b1=-24l-12l^2+r(16m+8m^2-16)
b2=12l^2+24l-12+r(24l^2-48lm+16m^2+16-16m).
```

The constant `12` is the discriminant at `k=0` for the target `y=-2`.
With the standard quadratic Sylvester determinant,

```text
Res_k(D4,D-2)
 =12(12 a2^2-b1 a1 a2+b2 a1^2)
 =4608 r(l-m)^2 F(l,m,r).
```

Thus, after imposing the degree-two conditions `r != 0` and `l != m`, the
displayed `F=0` is exactly the projective resultant condition for the two
target discriminants.  It is a candidate condition for a member with branch
values `4,-2`; one must still reconstruct the parameter and check that the
resulting `Pk,Qk` are coprime and define a nonconstant degree-two map.

The word “affine” matters.  Homogenizing `k=[K:L]` gives

```text
D4h  = K(a1 L+a2 K),
D-2h = 12L^2+b1 K L+b2 K^2.
```

When `a2=b2=0`, the resultant can vanish at the projective point `k=infinity`
even though there is no finite affine solution.  Two rational examples are:

```text
(l,m,r)=(0,1,3/4):    F=0,  D4=36k,  D-2=12+6k.
```

There is no finite common parameter, and the point at infinity has
`P_infinity=0`, hence is a constant map and must be discarded.

```text
(l,m,r)=(-2,-4,1/4):  F=0,  D4=48k,  D-2=12+12k.
```

Again there is no finite common parameter, but the point at infinity is a
valid degree-two map and does have branch values `4,-2`.  Consequently an
affine statement must exclude or separately handle `a2=b2=0`; a projective
statement must retain that locus and perform the degree/coprimality check.

For rational `l,m,r`, when `a2 != 0` the finite candidate is the rational
value `k=-a1/a2`.  If `a2=0` and `b2 != 0`, `F=0` forces `a1=0`, so `D4` is
identically zero and the roots of `D-2` still need to be reconstructed (and
need not be rational without a square check).

The representation scale explains the parameter in the `l=2` component.
If `t` is the coefficient of the unscaled pair `(A,A-rB)` relative to
`(P0,Q0)`, then `t=kc/(1-k)`. On the nondegenerate `l=2` component,

```text
r=9/[4(m^2-m+1)],  t=1/3,
k=1/(3c+1).
```

The displayed finite formula requires `c != 0` and `3c+1 != 0`; for
`c=-1/3` the same projective member is the point `k=infinity`.

Finally, substituting `l=2` and the displayed `r` into the common-argument
quartic gives the exact discriminant

```text
-729(m-2)^4(53m^2-50m+50) / [2(m^2-m+1)^3].
```

It is negative for every real `m != 2`, since `m^2-m+1>0` and
`53m^2-50m+50>0`.  At `m=2` the quartic has a repeated `(x-2)^2` factor;
this is also the degenerate `l=m` case.  Therefore this component cannot
give four distinct rational common arguments.
