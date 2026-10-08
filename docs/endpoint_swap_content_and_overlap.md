# Exact endpoint-swap contents and the remaining overlap cost

The endpoint involution has an exact coordinate-content formula for arbitrary
primitive rational endpoints. It gives a simultaneous denominator-clearing
parameter. A separate sufficient radius bound retains the full Gaussian
overlap of all source denominators. Neither statement bounds that overlap
uniformly or supplies a suitable map for arbitrary large endpoint clusters.

## 1. General endpoint and exact ordinary content

Let `x,y` be positive coprime integers. The endpoints have half-angle rows
`(1,0)` and `(x,y)`; a row `(u,v)` represents the circle phase
`(u+iv)/(u-iv)`. For a positive integer `D`, set

```text
M_D = [[xy,D-x^2],[y^2,-xy]].
```

Its square is `y^2 D` times the identity. On the half-angle interval
`0<=r<=y/x`, its action is

```text
r -> y(y-xr)/(xy+(D-x^2)r).
```

The denominator is positive: its endpoint values are `xy` and `Dy/x`.
The derivative is negative, and the endpoints are exchanged. Thus the
image retains the exact angular span `2 arctan(y/x)`. The map is conformal
exactly when `D=x^2+y^2`; other positive parameters are nonconformal.

For any primitive integer row `(u,v)`, including either endpoint, put

```text
d=yu-xv,       c=gcd(d,v)=gcd(y,v)>0,
d'=d/c,       v'=v/c,       h=gcd(d',D)>0.
```

Here gcds are positive, allowing negative or zero arguments. Its image is
`(xd+Dv,yd)`, with exact ordinary coordinate content

```text
content(M_D(u,v))
  = c h gcd(x(d'/h)+(D/h)v', y).                       (1)
```

Indeed `gcd(d',v')=1`. Write `d'=h a,D=h b`, so `gcd(a,b)=1`.
After removing `ch`, the image is `(xa+bv',ya)`, and
`gcd(a,xa+bv')=1`. Its remaining gcd is therefore exactly
`gcd(xa+bv',y)`. This argument includes `d'=0`: then `a=0,b=1`
and `v'=+1` or `-1`. It also includes the anchor `v=0`.

In particular, when `y=1`, equation (1) becomes

```text
content(M_D(u,v)) = gcd(u-xv,D).                       (2)
```

If `D` is divisible by every positive `d'` of the interior rows, their
primitive image rows are

```text
( (x+(D/d')v')/j, y/j ),
j=gcd(x+(D/d')v',y).
```

All finite image denominators consequently divide `y`. This is ordinary
coordinate normalization, not yet the full Gaussian circle normalization.

## 2. The exact extra multiplier for integer cotangents

Suppose `L>0`, `A=min_i X_i>0`, and distinct integer finite coordinates
`X_i` satisfy `(X_i X_j+L^2)/(X_i-X_j)` integral for every pair.
Put

```text
g=gcd(A,L),       x=A/g,       y=L/g,       T=x^2+y^2,
g_i=gcd(X_i,L),   H_i=(X_i/g_i,L/g_i).
```

For an interior row `X_i>A`, the quantity from Section 1 is exactly

```text
d_i'=(X_i-A)/gcd(g,g_i).
```

This follows from
`gcd(L/g,L/g_i)=L gcd(g,g_i)/(g g_i)`. Pair integrality with `A`
implies `X_i-A | A^2+L^2=g^2 T`. Define

```text
E=lcm_i [d_i'/gcd(d_i',T)].
```

Empty lcms are one. Primewise,

```text
v_p(E)=max_i max(0,v_p(d_i')-v_p(T)),
E | g^2,
T E = lcm(T,d_1',d_2',...).                            (3)
```

Thus `D=TE` is the least multiple of `T` clearing all these transverse
determinants. If it is conformal, `D=2TE` is a nonconformal clearing
parameter. No radius bound follows merely from this choice. All `d_i'`,
`T` and `E` are invariant under common scaling of `A,X_i,L`; the upper
bound `g^2` can become looser under scaling.

## 3. Full Gaussian radius normalization

For Gaussian half-angle representatives `H_0=1,H_1,...,H_s`, reduce each
phase fraction `H_i/bar(H_i)` in `Z[i]`. Let `A_i` be its reduced
denominator and `Lambda=lcm_G(A_1,...,A_s)`, defined up to a unit.
The least primitive squared circle radius is exactly

```text
N=Norm(Lambda).                                       (4)
```

Clearing these denominators constructs integral equal-norm rows, including
the anchor `Lambda`. Their Gaussian gcd is a unit: at every prime dividing
`Lambda`, a denominator attaining its maximal valuation has numerator a
unit at that prime, and its cleared row has valuation zero. Conversely,
every integral anchor clearing all phase fractions is divisible by
`Lambda`. This proves (4), with no replacement by the ordinary lcm of
the anchor-denominator norms.

For any fixed Gaussian denominator `B` and further denominators `A_i`,
put `Q_i=A_i/gcd_G(B,A_i)`. Then the exact identity is

```text
lcm_G(B,A_1,...,A_r) = B lcm_G(Q_1,...,Q_r)             (5)
```

up to a unit. In particular `Q_i` need not be coprime to `B`.
At a Gaussian prime with valuations `b,a_i`, its valuation is
`max(a_i-b,0)`; adding `b` to their maximum proves (5), including excess
powers over primes already dividing `B`.

## 4. A source-aware sufficient overlap bound

Here impose the following explicit hypotheses:

```text
y=1, x>0 even, T=x^2+1, B_0=x+i;
H_i=u_i+i v_i=d_i+v_i B_0,  i=1,...,r;
v_i>0, d_i=u_i-xv_i>0, gcd(d_i,v_i)=1;
k>=1 an integer, gcd(k,T)=1, D=kT.
```

The full source includes `1,B_0` and all the distinct interior phases
represented by the `H_i`. The case `k=1` is conformal; `k>1` is not.
Set `C_i=d_i+k v_i B_0`. The transformed representative satisfies

```text
M_D H_i = d_i B_0+k T v_i = B_0 bar(C_i).
```

After dividing all image phases by the phase of `B_0`, the image endpoint
is `bar(B_0)/B_0`, and interior phases are `bar(C_i)/C_i`.
Hence all target denominators come from `B_0` and the conjugate-primitive
reductions of the `C_i`. This reanchoring is at an actual image point.

Each `H_i` has relatively prime ordinary coordinates; `C_i` can have
ordinary content `gcd(d_i,k)`. Write `epsilon_i,epsilon_i'` for the norms
of their respective Gaussian gcds with their conjugates. Thus `epsilon_i`
is one or two, while `epsilon_i'` can be larger. Conjugate all source
denominators if necessary, so their anchor denominator is `B_0`.
Let `A_i` denote these reduced source denominators, and `A_i'` the reduced
target denominators. They have norms

```text
Norm(A_i)=Norm(H_i)/epsilon_i,
Norm(A_i')=Norm(C_i)/epsilon_i'.
```

The exact anchor gcd norms are

```text
Norm(gcd_G(B_0,A_i))=Norm(gcd_G(B_0,A_i'))=gcd(d_i,T).         (6)
```

To see this, `B_0` is coprime to its conjugate since `x` is even. At each
prime `pi|B_0` above `p`, let `e=v_pi(B_0)=v_p(T)` and `t=v_p(d_i)`.
If `t>0`, `v_i` and `k` are units at `p`. The valuation of `d_i+v_iB_0`
is `t` when `t<e`, is `e` when `t>e`, and is at least `e` when `t=e`;
its gcd with `B_0` therefore has valuation `min(t,e)`. The same holds
with `k v_i`. Their conjugates are units at `pi`, so reduction by the
full conjugate gcd cannot remove these factors, even when `C_i` has
ordinary content elsewhere. For `t=0` both expressions are units at `pi`.
This proves (6) with all common contents and multiplicities retained.

Put

```text
Q_i=A_i/gcd_G(B_0,A_i),
Q_i'=A_i'/gcd_G(B_0,A_i'),
Omega= product_i Norm(Q_i) / Norm(lcm_G(Q_1,...,Q_r)).
```

The positive integer `Omega` retains every shared Gaussian prime power.
Equations (4)--(6) give the exact source and target formulas

```text
N  = T Norm(lcm_G(Q_1,...,Q_r)),
N' = T Norm(lcm_G(Q_1',...,Q_r')).                     (7)
```

Since all `d_i,v_i,x` are positive and `k>=1`,

```text
Norm(C_i)=d_i^2+2kx d_i v_i+k^2 T v_i^2
          <= k^2 Norm(H_i).
```

Thus `Norm(Q_i')<=2k^2 Norm(Q_i)`. Taking a product only on the target
side of (7) proves the sufficient estimate

```text
N' <= (2k^2)^r Omega N.                               (8)
```

If `k` is odd, `H_i,C_i` have identical coordinate parities, so
`epsilon_i' >= epsilon_i`, improving the factor to `k^(2r)`. Because angular
width is preserved exactly, the corresponding normalized widths satisfy
`C' <= [(2k^2)^r Omega]^(1/4) C`, with the same odd-`k` improvement.

This is a sufficient inequality, but the product budget has the obstruction
in Section 5 and does not provide a viable uniform certificate for five or
more points at shrinking angular width. No condition `d_i|T` is needed
in this overlap estimate; it is separate from the denominator-clearing
construction in Section 2. The radius formulas keep the full Gaussian
lcms; no disjoint-support or squarefree assumption was made.

## 5. The product certificate already diverges for five points

The estimate in (8) discards target overlap. That loss is decisive, even
when the original source width tends to zero. Equation (6) gives

```text
Norm(Q_i)=Norm(H_i)/(epsilon_i gcd(d_i,T))
         >= Norm(H_i)/(2d_i)
         = d_i/2+xv_i+T v_i^2/(2d_i)
         >= (x+sqrt(T))v_i >= 2x.
```

The penultimate inequality is the arithmetic--geometric mean inequality;
the last uses `v_i>=1` and `T=x^2+1`. Write `Delta=2 arctan(1/x)`.
Since `Delta>=1/x` and `T>=x^2`, the exact source product budget obeys

```text
Delta^4 Omega N = Delta^4 T product_i Norm(Q_i)
               >= 2^r x^(r-2).                       (9)
```

For every fixed `r>=3` (at least five points including the endpoints),
this diverges as `x` tends to infinity. Multiplying it by the additional
factor in (8) only increases it. Thus this particular upper-bound
certificate cannot prove a fixed target endpoint constant along such a
sequence, regardless of the source's small endpoint constant. Equation
(9) is not a lower bound on the actual target radius: target denominator
overlap may still save the map. The four-point case `r=2` is the critical
boundary not excluded by this calculation, consistent with the Pell
repair. A useful general estimate must retain more target overlap than
the product step in (8).

There is also a necessary quantitative target-overlap bound. Define

```text
gamma_i=gcd(d_i,k),       w_i=(k/gamma_i)v_i,
Omega'_k=product_i Norm(Q_i')/Norm(lcm_G(Q_1',...,Q_r')).
```

The ordinary primitive representative of `C_i` is
`C_i/gamma_i=(d_i/gamma_i)+w_i B_0`. Its two positive coefficients
are coprime. Moreover `gamma_i` is coprime to `T`, so
`gcd(d_i/gamma_i,T)=gcd(d_i,T)`. Applying the preceding norm argument
to this primitive representative gives

```text
Norm(Q_i') >= (x+sqrt(T)) w_i >= 2x w_i.
```

Consequently, with the actual target normalized width `C'_k`,

```text
Omega'_k (C'_k)^4
  = Delta^4 T product_i Norm(Q_i')
  >= 2^r x^(r-2) product_i w_i.                        (10)
```

If `C'_k<=K`, its denominator overlap must therefore satisfy
`Omega'_k >= 2^r x^(r-2) product_i w_i / K^4`. Unlike (9), this
is a necessary condition on a successful target, rather than a limitation
of an upper-bound certificate. It does not assert that such overlap is
impossible. It retains every common Gaussian prime power and every
ordinary image content, including those outside `T`.

The [simultaneous alignment theorem](endpoint_swap_simultaneous_alignment_capacity.md)
uses Gaussian resultants to bound the product of these overlaps over
different fixed integer parameters. It limits simultaneous preservation
without excluding one successful parameter.

## Verification and scope

The companion [checker](check_endpoint_swap_content_and_overlap.py)
checks signed general contents, the integer-cotangent multiplier,
ramified reductions, shared powers remaining in `Q_i`, and the exact
all-row source/target radii by independent primitive tuple construction
and an ordinary lcm using every edge denominator. It passes 668,160
general content cases, 2,182 integer-cotangent cliques, and 199 all-row
overlap cases. A separate agent audited equations (1)--(8) and replayed
the checker without finding a correction. Root removed the unnecessary
condition `d_i|T` and independently audited the generalization. Another
Astra instance derived the product-budget obstruction and necessary
target-overlap bound; root independently checked both and replayed the
completed checker.
The finite checks supplement the proofs. These lemmas do not improve the
general circle-point bound or prove uniform endpoint-map existence.
