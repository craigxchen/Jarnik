# Isometric virtual anchors preserve radius but can acquire the missing content

This is a cost audit of the proposed extension of the
[actual-anchor shear obstruction](parabolic_shear_anchor_content_obstruction.md).
It gives an exact rational normalization, an exact virtual-anchor conductor
formula, and a Pell family saturating the content cost. It does not exclude
general critical maps or construct configurations of unbounded cardinality.

## 1. Exact rational normalization

Let `M=((a,b),(0,d))`, `delta=ad>0`, and

```text
D=b^2+(a-d)^2,       K=D(D+4delta).
```

A nonzero real vector `H=(x,y)` is isometric when
`Norm(MH)=delta Norm(H)`. Its binary quadratic equation is

```text
a(a-d)x^2+2abxy+[b^2+d(d-a)]y^2=0.
```

For a nonconformal map it has rational projective roots precisely when
`delta D` is a rational square, hence an integer square `s^2`.
Put `w=MH`, and use complex multiplication matrices `R_H` and `L_w`
for multiplication by `H` and `conjugate(w)`. Then exactly

```text
L_w M R_H = Norm(H) ((delta, +/-s),(0,delta)).          (1)
```

The upper right entry follows either by direct multiplication or from
`T/delta=2+(off_diagonal/diagonal)^2`. Thus, after primitive scalar
reduction, the shear is `((A,+/-B),(0,A))`, where

```text
A=delta/gcd(delta,s),       B=s/gcd(delta,s).
```

There is no surviving factor `Norm(H)` in its matrix height. The existing
source and image relative phase configurations, without adjoining the
virtual point, retain their exact least circle norms: all pair ratios
are unchanged by a common rotation. Here a realization allows a common
complex multiplier; the reference phase need not itself be an integer
point unless it belongs to the tuple. This does not make the virtual
direction an actual source point. Adjoining it has the separate cost
computed next.

When `a!=d`, the two root representatives can be taken as

```text
H_+ = (+s-ab, a(a-d)),       H_- = (-s-ab, a(a-d)),
```

then reduced by their separate ordinary contents. Their unreduced norm
product is `a^2(a-d)^2 K`. This identity includes both roots and their
individual normalization costs.

## 2. Exact radius and content of an adjoined rational anchor

Start from the odd primitive source setup of the shear note, anchored at
an actual point `j`. At each odd split prime let the signed source phase
interval be `J_p=[L_p,U_p]`, of width `e_p=v_p(N)`. The actual anchor
has signed coordinate zero, so `L_p<=0<=U_p`. For a primitive Gaussian
integer half-angle vector `H`, put

```text
r_p=v_pi(H)-v_bar(pi)(H),
E_H=product p^dist(r_p,J_p),
G_H=product p^max(0,min(r_p-L_p,U_p-r_p)).
```

The exact least squared radius after adjoining that direction is

```text
N_plus=N E_H.                                        (2)
```

In a primitive physical realization of this augmented configuration,
the ordinary coordinate gcd of the virtual point is exactly `G_H`.
Indeed its two split valuations are
`r_p-min(L_p,r_p)` and `max(U_p,r_p)-r_p`; their minimum is the exponent
displayed in `G_H`. An even half-angle norm contributes no factor two
to the least circle norm: its common `1+i` cancels in `H/conjugate(H)`.
Pair half-angle norms can retain a parity factor; the least circle
conductor does not. Consequently `gcd(E_H,G_H)=1`.

Under a rational isometric reanchoring, the old clipped conductor is
unchanged: source and coefficient intervals translate by the same `r_p`.
Applying the shear interval proof, either directly to the translated
interval or to the augmented tuple, yields

```text
Nclip | B(B^2+4A^2) G_H.                             (3)
```

The direct proof does not require the virtual point to be adjoined
geometrically. At a prime dividing `B`, the translated coefficient interval
is `[-t,t]`; its intersection with the translated source interval has
width at most `t+max(0,min(r_p-L_p,U_p-r_p))`.

The general comparison with the actual anchor content is only

```text
G_H | g_j oddpart(Norm(H)),                          (4)
```

by the one-Lipschitz dependence of interior depth on the signed coordinate.
Conductor sharing can make `E_H=1` while making `G_H` large. Thus removing
the radius cost does not remove the content cost.

Both isometric roots can incur the same cost. For an odd prime
`p|s`, `p` not dividing `2ad(a-d)`, set `t=v_p(s)`. The original
coefficient interval is either `[0,2t]` or `[-2t,0]`, and **both**
isometric directions have signed coordinate its midpoint, `+t` or `-t`.
To verify both orientations, write `V=(a-d)+ib` and
`H_+ = s+iaV`, `H_-=-s+iaV`. Exactly one orientation of `V` has
valuation `2t`; each `H` has valuation `t` in that orientation and zero
in the other. Their ordinary primitive reductions are units at `p`.

## 3. An exact Pell family saturates the loss

Take positive odd solutions `b^2-2c^2=-1`, with `c>=5`, and set

```text
M=((1,b),(0,2)),       N=c^2(c^2+4),
z_0=-(c^2+1)-ib,      z_1=-(c^2+1)+ib.
```

Both original physical points have ordinary content one and norm `N`.
Their Gaussian gcd is a unit, so their exact least circle norm is `N`.
For the content claim, any common divisor of `b` and `c^2+1`, or of
`b` and `c^2-2`, divides three. The Pell equation modulo three forces
both `b` and `c` to be units there. Parity then also excludes two from
the target coordinate gcd. A conjugate pair of coordinate-primitive
Gaussian integers of odd norm has Gaussian gcd one.
The relative half-angle rows are `1` and `conjugate(z_0)`.
Their images have primitive second row `(c^2-2,2b)`, also of norm `N`;
the image least circle norm is exactly `N`. Explicitly, the common target
reference is `w_0=c^2-2-2ib`. Thus the source phase is `z_1/z_0`, and
the image phase is `conjugate(w_0)/w_0`. Physical virtual points are
`z_0 H/conjugate(H)` in the source and
`w_0 (MH)/conjugate(MH)` in the target; they are not scalar multiples
of their half-angle rows.

Here `s=2c`, and the two isometric rows are

```text
H_+=(-b+2c,-1),       H_-=(-b-2c,-1).
```

Their physical source points and physical image points are respectively

```text
v_+=-c^2+2ic,         v_-=-c^2-2ic,
w_+= c^2+2ic,         w_-= c^2-2ic.                  (5)
```

In the image realization the original two points are
`(c^2-2)+/-2ib`. All four source points and all four image points
are integral on their existing circles. Adding either or both virtual
anchors leaves both least norms exactly `N`. Each virtual point has
ordinary content `c`, on both sides. The primitive parabolic matrices
in (1) have `A=1, |B|=c`.

At every source prime the two original signed phase values cover the
coefficient interval. To check this without a prime-orientation shortcut,
use `U=3-ib`, `V=-1+ib`, and `z_0=U conjugate(V)/2`.
Their norm product is `4N`; at odd primes neither Gaussian coefficient
has both orientations, and their odd norm parts are coprime, since
`gcd(c^2,c^2+4)=1`. The source extrema are exactly the two coefficient
endpoints. Therefore

```text
Nclip=N=c^2(c^2+4)=B(B^2+4) G_H.                    (6)
```

This saturates (3) with `g_0=g_1=1`, `G_H=c~N^(1/4)`, and zero
additional radius cost. The norms of both `H_+` and `H_-` are of order
`c^2`, so virtual height itself is of order `N^(1/2)`.

The example also keeps bounded endpoint geometry. Both four-point tuples
lie on arcs of length at most `8 N^(1/4)` for `c>=5`: the source angular
width is at most `2 arctan(2/c)`, and the target width is at most
`2 arctan(2b/(c^2-2))`. Thus the cost barrier is not caused by a distant
virtual direction. This is a four-point example, not an obstruction to a
theorem that uses arbitrarily many points or requires smaller constants.
The uniform arc certificate uses `N^(1/4)<=6c/5`, `b<=3c/2`,
`c^2-2>=9c^2/10` for `c>=5`, and `arctan(t)<=t`.
The original source pair's normalized arc constant tends to `2sqrt(2)`,
and the original image pair's tends to `4sqrt(2)`. In particular, this
family does not meet the `C,C'<=2` hypotheses of the actual-anchor
obstruction. A refinement using those smaller constants and many points
remains open.

## 4. Irrational isometric directions do not give a rational shear

If `delta D` is nonsquare, (1) is instead over
`L=Q(sqrt(delta D))`. The two isometric rows and the two signs of the
shear are Galois conjugates. A rational equal-diagonal representative
under conformal input/output changes is impossible: its rational shear
ratio would have square `D/delta`, contradicting nonsquareness.
Returning to rational coordinates at an actual source point therefore
cannot retain that equal-diagonal form.

The Gaussian integer clipping proof cannot simply be applied to these
algebraic entries. It would need valuations in `L(i)`, with ramification
and residue degrees, and an algebraic replacement for `G_H`. In particular,
the two conjugates cannot be counted as independent rational anchors.
Even the formal numerator has the exact field norm

```text
|Norm_(L/Q)(s(s^2+4delta^2))|=delta^3 D(D+4delta)^2.
```

This records the conjugation cost; it is not a new rational divisibility.
No algebraic-integer clipping theorem or bound on its anchor content is
asserted here. The rational Pell family already shows that rationality,
exact conductor sharing, and both isometric anchors together do not force
the small content needed by the actual-anchor obstruction.

The accompanying [exact checker](check_isometric_virtual_anchor.py) verifies
the normalization, both roots, least norms, contents, full clipped widths,
and the sharp Pell identities.
