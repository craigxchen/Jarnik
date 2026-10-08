# Composition and height in the rational endpoint pencil

This note records an exact group law and a source-independent inequality
between the primitive radii of two configurations in one rational endpoint
orbit. The four-row pairing below is the half-angle form of the existing
[four-seed secant involution](fixed_secant_conductor_bridge.md), not a new
route to full-row integrality. The height inequality applies to **all retained
rows**: every radius is the least
primitive squared radius computed from the full Gaussian denominator lcm,
equivalently from the ordinary lcm of every pair-edge norm. The proof of the
inequality uses one edge as a divisor of that full radius; it never replaces
the full lcm by a product of row norms.

Let `x>0` be an integer, `T=x^2+1`, and `B=x+i`. Take the endpoints `1,B`
and at least one distinct interior primitive half-angle row

```text
H_i=d_i+v_i B,    d_i,v_i>0,    gcd(d_i,v_i)=1.
```

The interior condition ensures that its phase lies strictly between the
endpoint phases. For coprime positive integers `p,q`, define

```text
G_(p/q) = [[q, x(p-q)], [0,p]],
G_(p/q) H_i = q d_i+p v_i B.
```

Both endpoints are fixed projectively, and the shortest containing angular
span `Delta=2 arctan(1/x)` is unchanged. Write `N_(p/q)` for the least
primitive squared radius of the **full** image tuple. The parameter `1` is
the identity and every other positive parameter is nonconformal.

## 1. Exact composition and pair reciprocity

The matrices multiply projectively according to

```text
G_mu G_lambda ~ G_(mu lambda),       G_lambda^(-1) ~ G_(1/lambda).
```

For the endpoint-swapping involutions

```text
M_(p/q)=[[qx,pT-qx^2],[q,-qx]],
M_1=[[x,1],[1,-x]],
```

one has the exact integral identity

```text
M_1 M_(p/q) = T G_(p/q).
```

The map `M_1` is a conformal circle reflection. Thus swap and fixing
parameters have the same full primitive target radius.

Put `w_i=d_i/v_i`. The fixing map sends `w_i` to `w_i/lambda`, while the
conformal reflection sends `w_i` to `T/w_i`. For any two interior labels
`i,j`, the positive rational parameter

```text
lambda_ij = w_i w_j/T
```

therefore maps the four-row subset `{1,B,H_i,H_j}` to its conformal
reflection, exchanging the two interior labels. Its **exact** four-row
primitive squared radius equals the source four-row radius. It is
nonconformal when `w_i w_j != T`.

This argument extends to the full tuple precisely when the multiset of
interior `w_i` is closed under `w -> lambda T/w`: then `G_lambda` maps the
full phase set to its conformal reflection and `N_lambda=N_1`. Pairwise
reciprocity alone gives no bound on the other rows' denominators. In
particular, composing subset symmetries does not produce a full-row
radius bound.

## 2. A proper height inequality for the full orbit

For every pair of positive rational parameters `lambda,mu`, write the
relative parameter `mu/lambda=P/Q` in lowest positive terms. Then

```text
             x(P+Q) <= 2 sqrt(N_lambda N_mu).          (1)
```

To prove this, regard the primitive full image at `lambda` as the source.
Its endpoints are still `1,B`, and any interior phase has a primitive
representative `H=d+vB` with `d,v>0` and `gcd(d,v)=1`. Let `N=N_lambda`.
The reduced Gaussian denominator norm of the edge from `1` to `H`
divides `N`. Its ordinary coordinate content is one, and possible
Gaussian ramification removes at most a factor two. Hence

```text
d+xv <= sqrt(2N).                                      (2)
```

The relative map sends this row to `C=Qd+PvB`. Its ordinary coordinate
content is exactly

```text
gamma=gcd(Qd,Pv)=gcd(P,d) gcd(Q,v) <= d v.            (3)
```

The equality follows prime by prime from `gcd(P,Q)=gcd(d,v)=1`. The
reduced Gaussian denominator norm of the target edge from `1` to `C`
divides `N_mu`. Again allowing the factor two at `1+i` gives

```text
sqrt(2N_mu) >= (Qd+Pxv)/gamma
              >= Q/v + Px/d
              >= x(P+Q)/sqrt(2N),
```

where the last step uses (2) separately for `d` and `xv`. This proves
(1), with every ordinary content and ramified cancellation included.

For `lambda=p/q` and `mu=r/s` in lowest terms, equation (1) reads

```text
x (rq+sp)/gcd(rq,sp) <= 2 sqrt(N_(p/q) N_(r/s)).      (4)
```

Taking `lambda=1` shows that, for a fixed source squared radius `N`, every
parameter whose full image has squared radius at most `B_0` satisfies

```text
p+q <= 2 sqrt(N B_0)/x.                               (5)
```

Thus the bounded-radius part of this rational orbit is finite, with at
most `O(N B_0/x^2)` parameter pairs. This is a source-dependent bound.
For endpoint-scale families with `N` and `B_0` both of order `x^4`, it
allows parameter heights of order `x^3` and does not provide a uniform
number of maps or a nonconformal image. Reversibility only says that a
chosen image can be mapped back; it cannot turn the identity's bounded
radius into a second bounded point of the orbit.

The same estimate has a rational-endpoint form. If `a,b>0` are coprime,
`B=a+ib`, and `H=d+vB` with coprime positive `d,v`, then

```text
a(P+Q) <= 2b^2 sqrt(N_lambda N_mu),                  (6)
```

or, with `x=a/b`, `x(P+Q)<=2b sqrt(N_lambda N_mu)`.
The ordinary coordinate content of `H` is at most `b`, so its edge norm
gives `d+av<=b sqrt(2N_lambda)`. The image `Qd+PvB` has content at most
`b gcd(Qd,Pv)<=bdv`; repeating the proof of (1) yields (6). This
includes every possible prime dividing `b` without an assumption that
the displayed `H` is coordinate-primitive.

There is an exact finite no-descent example within this pencil. Take
`x=1`, endpoints `1,1+i`, and three interior rows
`2+i,3+i,7+i`. Their full primitive squared radius is `125`.
For a parameter `p/q`, the first interior row becomes `(p+q)+ip`,
whose ordinary content is one. If the full target squared radius were at
most `125`, its anchor-edge norm would imply
`(p+q)^2<=250`, hence `p+q<=15`. An all-edge lcm check of the 71
coprime positive pairs in this finite range finds `N_(p/q)<=125` only
at `p=q=1`. Thus identity is a strict global radius minimum among
positive rational parameters for this five-row tuple. This does not
exclude a nonconformal image at a larger fixed endpoint constant;
`p/q=3` has squared radius `425`.

The [checker](check_endpoint_swap_composition_orbit_height.py) verifies
the exact matrix and pair symmetries, ordinary contents, primitive-radius
identities, both height inequalities, and the finite no-descent fixture.
