# A height floor for adaptive endpoint alignments on the cubic five-point family

The [cubic five-point endpoint family](endpoint_pair_selected_five_point_obstruction.md)
defeats every fixed nonidentity positive rational alignment with its two
physical endpoints retained. Here the rational parameter is allowed to
vary with the family parameter. A bound on the **true primitive radius**
of a three-row target subset forces its reduced parameter height to grow.
This restricts adaptive choices on this one actual endpoint family; it
does not rule them out and does not prove a general circle-point bound.

Let `t` be a positive multiple of `600`. Use the source rows `1,B,H_1,H_2,H_3`
and the exact endpoint angle `Delta` from the cited family. Its least
primitive squared source radius satisfies `N~t^12/720^2`, and its
normalized width tends to `2 sqrt(5)`. Write

```text
B=P+iQ,       P=t^4+25t^2-36,       Q=60t,
d=20(t^2+4)(t^2+9)(t^2+36),       v=40(t^2+18).
```

For coprime positive integers `p,q` with `p!=q`, put `H=p+q` and use the
endpoint-fixing rational alignment `lambda=p/q`. Its third interior
target row is `C=qd+pvB`; the full target retains all five labels.
Let `N'_(p/q)` denote its least primitive squared radius.

**Uniform adaptive-parameter bound.** With the absolute integer

```text
K0=40*45360*60*25920=2821754880000,
```

every such `t,p,q` satisfies

```text
N'_(p/q) >= t^14 H^2/[62208 K0^2 p^2 q^3 (p-q)^2]
          >= t^14/[62208 K0^2 H^3 (p-q)^2]
          >= t^14/[62208 K0^2 H^5].                  (1)
```

The first bound retains the relative balance of the numerator,
denominator, and their difference. It is derived from the Gaussian
denominator lcm of the actual subset `1,B,C`, whose primitive squared
radius divides the full target radius. No row-norm product is used as a
substitute for that lcm.

## Exact content and overlap bounds

Set

```text
U=(t+2i)(t+3i)(t-6i),       S=t+i.
```

Then, exactly in `Z[i]`,

```text
B=US,       d=20 U bar(U),
C=U E,      E=20q bar(U)+pvS.                         (2)
```

For an integral half-angle row `Z`, let `D_Z` be the denominator of its
reduced Gaussian phase fraction `Z/bar(Z)`. If `g_Z` is its ordinary
coordinate content, then

```text
Norm(D_Z)=Norm(Z)/(g_Z^2 epsilon_Z),
epsilon_Z in {1,2}.                                 (3)
```

The factor two is the complete possible cancellation at `1+i` after
ordinary primitivization. On multiples of `600`, the ordinary content
of `B` is exactly `36`: `P==-36 mod 60t` and `36|60t`. The coordinates
of `B/36` are odd and even, so `epsilon_B=1`. Since `Norm(B)>=t^8`,

```text
Norm(D_B)=Norm(B)/36^2 >=t^8/36^2.                  (4)
```

Let `gamma=gcd(Re C,Im C)` and `A=Re C=qd+pvP`. As `Im C=pvQ`,
the ordinary gcd inequalities give

```text
gamma<=p gcd(A,vQ)<=p gcd(A,v) gcd(A,Q).              (5)
```

Now `gcd(A,v)=gcd(qd,v)<=q gcd(d,v)`. The polynomial identity
`d==45360 mod (t^2+18)` yields
`gcd(d,v)<=40*45360`. Also

```text
A==25920(q-p) mod t,
gcd(A,Q)=gcd(A,60t)<=60*25920*|q-p|.                 (6)
```

The last bound uses `p!=q`. Therefore

```text
gamma<=K0 p q |p-q|.                                 (7)
```

For `t>=600`, `P>=t^4`, `d>=20t^6`, and `v>=40t^2`. Thus
`A>=20(q+2p)t^6>=20H t^6`, so (3) gives

```text
Norm(D_C)>=Norm(C)/(2gamma^2)
          >=200 H^2 t^12/gamma^2.                    (8)
```

The remaining Gaussian overlap is controlled without assuming the
two reduced denominators coprime. Polynomial evaluation at `t=-i`
gives `bar(U)==-60i mod S`; hence from (2),

```text
gcd_G(B,C)=U gcd_G(S,E),
gcd_G(S,E) | 1200q.                                      (9)
```

Because `t` is even, `S=t+i` has odd norm and is coprime to its
conjugate. At each rational prime it uses at most one split Gaussian
orientation, while inert and ramified primes do not occur. Therefore,
for every positive rational integer `m`,

```text
Norm(gcd_G(S,m))=gcd(Norm(S),m)<=m.                    (10)
```

As `Norm(U)<=8t^6`, (9)--(10) yield

```text
Norm(gcd_G(B,C))<=9600q t^6.                         (11)
```

The reduced denominators `D_B,D_C` divide `bar(B),bar(C)`. Thus their
Gaussian gcd has norm at most the right side of (11), and their exact
lcm norm obeys

```text
N'_(p/q)>=Norm(lcm_G(D_B,D_C))
          =Norm(D_B)Norm(D_C)/Norm(gcd_G(D_B,D_C))
          >=t^14 H^2/(62208 gamma^2 q).              (12)
```

Substituting (7) proves the first inequality in (1); the others use
`p,q,|p-q|<=H`. This keeps all ordinary contents, split-prime
orientations, Gaussian gcd multiplicities, and the possible ramified
factor of `C`.

## Consequence for bounded endpoint images

For `t>=600`, one has `P<=2t^4` and `Q/P<=1`. Since
`2 arctan(u)>=u` on `[0,1]`, the exact endpoint span satisfies

```text
Delta=2 arctan(Q/P)>=30/t^3.                         (13)
```

If the full target normalized width obeys
`Delta (N'_(p/q))^(1/4)<=C'` for a fixed `C'>0`, (1) and (13)
force

```text
H^3(p-q)^2 >= [30^4/(62208 K0^2 (C')^4)] t^2.    (14)
```

In particular `H>=c_(C') t^(2/5)`, equivalently
`H>=c_(C') N^(1/30)`. If `|p-q|` is bounded along a sequence, then
`H>=c_(C',|p-q|) t^(2/3)`, equivalently
`H>=c_(C',|p-q|) N^(1/18)`. The identity alignment `p=q=1` is
excluded from (1)--(14): its ordinary content has a growing factor
and its target is the original bounded-endpoint tuple.

The [exact checker](check_endpoint_swap_adaptive_height_five_point.py)
compares independent Gaussian-lcm and ordinary all-edge primitive
radii while testing (2)--(12) on fixed and varying rational parameters.
The displayed integer arguments prove the bounds for all allowed
parameters and specializations; the finite checks supplement them.

## Scope of the all-row overlap calculation

For any interior row `C_h=qd_h+pv_h B`, the ordinary content has an
exact parameter factorization. Put

```text
g=gcd(d_h,v_h),       d0=d_h/g,       v0=v_h/g,
a=gcd(p,d0),         b=gcd(q,v0),
X=(q/b)(d0/a)+(p/a)(v0/b)P.
```

Then

```text
gcd(Re C_h,Im C_h)=gab gcd(X,Q).                     (15)
```

Indeed the two integers `(q/b)(d0/a)` and `(p/a)(v0/b)` are coprime.
After removing `gab`, the remaining content is therefore
`gcd(X,(p/a)(v0/b)Q)=gcd(X,Q)`.

For two interiors, write `W_ij=v_j d_i-v_i d_j`. Linear combinations
of their raw rows show that their Gaussian gcd divides both `qW_ij`
and `pBW_ij`. Since `gcd(p,q)=1`,

```text
gcd_G(C_i,C_j) | W_ij gcd_G(B,q).                    (16)
```

Including endpoint `B` in the full denominator lcm does **not** in
general recover the overlap factor `gcd_G(B,q)` as an extra gain.
Suppose an odd split prime has `v_pi(B)=e>0`, `v_bar(pi)(B)=0`, and
`v_p0(q)>e` at the underlying rational prime `p0`. The parameter
numerator `p` is a unit there. Every `C_h` then has valuations `e,0`
at these two orientations if `v_h` is also a unit at `p0`. Its
reduced phase denominator and that of
`B` all have the same orientation with exponent `e`; their lcm still
has exponent `e`, not twice that exponent.

This happens inside the present family: take `t=600`, parameter
`p/q=170/169`, and `pi=3-2i` above `13`. For `B` and both normalized
rows `qD_i+pB`, where
`D_1=4(t^2+1)(t^2+9)` and `D_2=(t^2+1)(t^2+36)`, the valuations
are exactly `1,0`. Their full three-row denominator lcm has exponent
one at `bar(pi)`. Moreover `13` does not divide
`D_1-D_2=3t^2(t^2+1)`. The checker verifies this exact prime-power
fixture in addition to (15)--(16) on its parameter cases. These
identities retain joint overlap information, but the resulting
all-row estimates have not improved (1) or excluded adaptive maps.
