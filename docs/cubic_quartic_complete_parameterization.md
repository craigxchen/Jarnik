# Complete quartic parameterization on the cyclic cubic norm cone

Let `K=Q(alpha)`, where

```text
alpha^3-3alpha+1=0,
w=2-alpha^2,       v=alpha^2+2alpha-2.
```

The trace form `q=e_2` satisfies

```text
Tr(w)=Tr(v)=0,     q(w)=-3,     q(v)=-9,
q(a+bw+cv)=3(a^2-b^2-3c^2).
```

Every quartic `U in K[t]` satisfying `q(U)=-3` and `U(i)=-i`
belongs to the two-parameter family described below. There are no such
polynomials of degree at most three. The parameters consist of a rational
rotation of the trace-zero basis and the rational parameter in the
[individual quartic construction](cubic_norm_descent_root_incidence.md).
An arbitrary rotation is a rational linear isometry of the trace form;
it need not preserve multiplication or splitting over `K(i)`.

## Coefficient reduction

Write `U=a+bw+cv`, with `a,b,c in Q[t]`, and put `s=t^2+1`.
The elements `1,w,v` remain independent over `Q(i)`, since the real
cubic field and `Q(i)` have coprime degrees. Thus `U(i)=-i` gives

```text
a=-t+s A,       b=s B,       c=s C,
deg(A),deg(B),deg(C) <= 2.
```

The equation `a^2+1=b^2+3c^2` becomes

```text
1-2tA+s A^2 = s(B^2+3C^2).                         (1)
```

Modulo `s`, this forces `A=-t/2 mod s`. Therefore, for rational `k`,

```text
A=-t/2+ks,
a=-(t^3+3t)/2+ks^2,
B^2+3C^2=(ks-t/2)^2+1-2kt.                         (2)
```

Set `beta_j=(B_j,C_j)`, where subscripts denote coefficients of
`t^j`, and use the inner product
`< (x,y),(x',y') > = xx'+3yy'`. Coefficient comparison in (2) gives

```text
<beta_2,beta_2> = k^2,
2<beta_2,beta_1> = -k,
<beta_1,beta_1>+2<beta_2,beta_0> = 2k^2+1/4,
2<beta_1,beta_0> = -3k,
<beta_0,beta_0> = k^2+1.                            (3)
```

If `k=0`, positivity of `x^2+3y^2` forces `beta_2=0`.
The remaining equations would give squared lengths `1/4` and `1`
for `beta_1,beta_0`, with inner product zero. Their Gram determinant
is both `1/4` and `3 det(beta_1,beta_0)^2`, impossible over `Q`.
This excludes all solutions of degree at most three as well.

## Normalizing the leading coefficient

Since `k!=0`, write `beta_2=k(p,q_0)`, where
`p^2+3q_0^2=1`. The matrix

```text
R = [[p, 3q_0], [-q_0, p]]
```

preserves `x^2+3y^2` and sends `beta_2` to `(k,0)`.
Apply it to the coefficient pair `(B,C)` and replace the basis by

```text
w'=p w+q_0 v,       v'=-3q_0 w+p v.
```

The resulting polynomial `U` is unchanged: the coordinate and basis
changes are inverse to each other. In these coordinates, suppressing
primes, the cubic coefficient in (2) gives `B_1=-1/2`.
Write `C_1=2ku`. If `u=0`, the quadratic coefficient would give
`B_0=k`, while the linear coefficient gives `B_0=3k`, a contradiction.
The remaining equations imply

```text
h=3u^2,
B_0=k(1-2h),
C_0=-(1+h)/(6u),
(h-1)(16k^2 h^2+h-1)=0.                            (4)
```

Because `3u^2=1` has no rational solution, the last equation yields
`16k^2 h^2=1-h`. Set `v_0=4kh`; then

```text
v_0^2+3u^2=1,       k=v_0/(12u^2),       u v_0 != 0.
```

Consequently every quartic has the form

```text
a=-(t^3+3t)/2+k(t^2+1)^2,
b=(t^2+1)[kt^2-t/2+k(1-6u^2)],
c=(t^2+1)[2ku t-(3u^2+1)/(6u)],
U=a+bw'+cv',
u=2r/(1+3r^2),       v_0=(1-3r^2)/(1+3r^2),
r in Q\{0},         p^2+3q_0^2=1.
```

Conversely, every indicated rational choice gives a quartic with the
required identities. With unrestricted rational `r`, the conic
parameterization misses `(u,v_0)=(0,-1)`; removing `r=0` also
removes `(0,1)`. Both points are excluded by `u!=0`, and every
valid point has inverse `r=u/(1+v_0)`. The denominator and `v_0`
do not vanish at any nonzero rational `r`.

The leading coefficient is `ell=k(1+w')`. It is nonzero because
`Tr(w')=0`, whereas `Tr(-1)=-3`. Since `q(ell)=0`, its inverse
`Y=1/ell` has trace zero. Thus `P=Y(U+i)` is monic with nonzero
constant imaginary part and constant imaginary norm. The element `Y`
is nonrational and therefore generates the cubic field. The full
family has both the conic parameter `r` and the independent rational
basis rotation; the fixed-basis family alone is not exhaustive.

## Rotation parameters at 17

Every rational rotation `(p,q_0)` is 17-adically integral. The form
`x^2+3y^2` is anisotropic over `F_17`, since `-3=14` is a nonsquare.
If `m=min(v_17(p),v_17(q_0))<0`, write
`(p,q_0)=17^m(P,Q)` with integral `P,Q`, at least one a unit.
Anisotropy makes `P^2+3Q^2` a unit, forcing
`v_17(p^2+3q_0^2)=2m<0`, contrary to its value one.
Hence every rational rotation reduces to one of the eighteen
norm-one points over `F_17`.

The [modulo-17 proof](cubic_norm_descent_mod17_nonsplitting.md)
therefore excludes every quartic in this classification for which
`v_17(r)=0`. Its fixed-basis theorem covers nonunit parameters only
for that basis. A basis isometry is not, in general, a field
automorphism, so it cannot transfer that theorem to all quartics.
