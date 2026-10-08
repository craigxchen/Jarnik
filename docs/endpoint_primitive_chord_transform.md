# The primitive-chord strip transform at the square-root endpoint

This note tests a route from the original integer circle equation that
does not use additive energy or bounded squared-difference
multiplicity.  The extreme chord of a short arc puts all of its lattice
points in a strip of constant Euclidean width.  Resolving that strip in
the primitive chord lattice gives an exact integral conic and a bound
that is strong when the endpoint chord has large coordinate gcd.

The transform is exact, but it does not prove a radius-independent
bound.  When the chord is primitive, its remaining count is of order
`sqrt(R)`, and the transformed conic is another sum-of-two-squares
problem at a worse normalized arc scale.

## 1. Direction normalization and the endpoint chord

Let the lattice points lie on

```text
x^2+y^2=n,       n=R^2 in Z,
```

in an arc of length at most `C sqrt(R)`.  As in
[the recent endpoint audit](endpoint_recent_literature_check.md),
multiplication by `1` or `1+i`, followed by coordinate sign changes,
is injective and puts the arc a fixed angular distance from both axes.
The new radius differs from `R` by at most `sqrt(2)`, so its two
coordinate intervals have length `O_C(sqrt(R))` and lower endpoints
comparable with `R`.  This normalization verifies the usual
simultaneous-square setup, but the chord transform below does not need
the separation from the axes.

Assume there are at least two points, and choose the first and last
points in cyclic order:

```text
P_0=(x_0,y_0),       P_1=P_0+e,
e=(u,v),             ell=|e|.
```

Let `L_0` be the arc length from `P_0` to `P_1`.  Then

```text
ell=2R sin(L_0/(2R))<=L_0<=C sqrt(R).                 (1)
```

For the finitely bounded range where `L_0/R` is not small, a bound may
depend directly on `C`; in the asymptotic range the chosen arc is the
minor arc.  Put

```text
g=gcd(|u|,|v|),       e_0=(u/g,v/g),
s=|e_0|^2=(u/g)^2+(v/g)^2,
G=gs=ell^2/g.                                         (2)
```

The vector `e_0` is primitive.  Equality of the endpoint norms gives

```text
2P_0 dot e+|e|^2=0,
P_0 dot e_0=-G/2.                                     (3)
```

In particular `G=gs` is even.

## 2. Exact dot-determinant coordinates

For another lattice point `P=P_0+h` on the arc define

```text
p=e_0 dot h,
q=det(e_0,h).
```

Both are integers.  Conversely,

```text
h=(1/s)(u_0p-v_0q, v_0p+u_0q),
e_0=(u_0,v_0),                                        (4)
```

so `(p,q)` comes from an integer `h` exactly when

```text
s divides u_0p-v_0q,
s divides v_0p+u_0q.                                  (5)
```

Thus the image is a sublattice of `Z^2` of index `s`; no integrality
condition is lost.

Set

```text
B=det(e_0,P_0).
```

The circle equation for `P_0+h`, together with (3)--(4), is exactly

```text
p^2+q^2-Gp+2Bq=0.                                    (6)
```

The endpoint coordinates are `(0,0)` and `(G,0)`.  Also

```text
(P_0 dot e_0)^2+det(e_0,P_0)^2=n s,
G^2+4B^2=4ns.                                         (7)
```

Completing squares in (6) gives the equivalent sum-of-two-squares
identity

```text
(G-2p)^2+4(q+B)^2=4ns.                                (8)
```

Equations (5), (6), and the box bounds below are a reversible
description of the original lattice points between the two endpoints.
They are also a backwards construction: choosing primitive `e_0`, an
integer `g` with `gs` even, and a lattice point `P_0` satisfying (3)
produces precisely this conic and sublattice.

## 3. Constant physical width becomes an arithmetic gcd bound

The maximum distance between a circular minor arc and its endpoint
chord is its sagitta

```text
W=R(1-cos(L_0/(2R)))
 <=L_0^2/(8R)
 <=C^2/8.                                             (9)
```

Since `|e_0|=sqrt(s)`, the perpendicular distance from `P` to the chord
line is `|q|/sqrt(s)`.  Every point of the minor arc is on the same side
of that line, and projection onto the chord lies between the endpoint
projections.  After changing the sign of `q` if needed,

```text
0<=p<=G,
0<=q<=W sqrt(s)
       <=(C^2/8)sqrt(s)
       <=(C^3/8g)sqrt(R).                             (10)
```

For a fixed integer `q`, equation (6) is quadratic in `p`, so it gives
at most two points even before imposing the sublattice congruences
(5).  Therefore the transform proves the explicit inequality

```text
M<=2(1+floor(W sqrt(s)))
 <=2+(C^2/4)sqrt(s)
 <=2+(C^3/4g)sqrt(R).                                 (11)
```

This uses convexity and the exact circle equation rather than an
energy count.  It is uniform whenever the extreme chord satisfies
`g >= c sqrt(R)` for some positive `c`, and it can be substantially
better than the trivial coordinate-interval count for an imprimitive
chord.

The exponent in the remaining case is explicit.  A primitive extreme
chord has `g=1`, `ell=O_C(R^(1/2))`,
`s=O_C(R)`, and (11) is only

```text
M=O_C(R^(1/2)).                                       (12)
```

There is no general lower bound for `g` in (3): when `e_0` is fixed,
the equation `P_0 dot e_0=-gs/2` is one integral linear condition, and
primitive choices of the chord direction are compatible with it.
This compatibility occurs at the endpoint scale, not only for tiny
fixed chords.  For `m>=1`, take

```text
e=(2m+1,2m-1),       g=1,       |e|^2=8m^2+2,
P_0=(-2m^2-m,2m^2+m+1).
```

Then `P_0 dot e=-(8m^2+2)/2`, so `P_0` and `P_0+e` have equal norm.
Their radius and chord length satisfy

```text
R asymp 2sqrt(2)m^2,       ell asymp 2sqrt(2)m,
ell/sqrt(R) -> 2^(3/4),                                  (13)
```

while the chord remains primitive.  This two-point family is not a
large-arc counterexample; it only rules out extracting a useful lower
bound for `g` from the equal-norm and endpoint-scale hypotheses.

## 4. The square-discriminant condition and its exponent loss

One might hope that requiring the quadratic (6) to have an integral
root for many values of `q` strengthens (11).  Its discriminant
condition is

```text
(G-2p)^2=4ns-4(q+B)^2,                                (14)
```

which is exactly (8).  Hence it is not a one-variable square condition
with a fixed small coefficient: it is another lattice-circle equation.

The affine transformation

```text
(p,q) |-> (r,t)=(G-2p,2(q+B))                         (15)
```

maps the points to a circle of radius

```text
R_*=2R sqrt(s)<=2C R^(3/2)/g.                         (16)
```

The dot-determinant map scales Euclidean lengths by `sqrt(s)`, and
(15) adds a factor two.  Thus the image arc has length

```text
L_*=2sqrt(s)L_0
    =2(ell/g)L_0
    <=2C^2 R/g.                                       (17)
```

In the difficult regime `ell` comparable with `sqrt(R)` and `g=1`,
these sizes are

```text
R_* asymp_C R^(3/2),       L_* asymp_C R.
```

Consequently

```text
L_*/sqrt(R_*) asymp_C R^(1/4),                        (18)
```

so the transformed arc is longer than the square-root endpoint by an
unbounded factor.  More generally the upper-scale version of this
factor is `O_C(R^(1/4)/sqrt(g))`.  When `g` is large enough to make that
factor bounded, (11) has already supplied a uniform count.  Iterating
the circle interpretation of (14) therefore gives no gain in the
uncontrolled primitive-chord case.

The sublattice conditions (5) have index `s`, but the rectangle in
(10) has horizontal length `G=gs` and vertical length
`O_C(sqrt(s))`.  Its area divided by the index is

```text
O_C(G sqrt(s)/s)=O_C(g sqrt(s))=O_C(ell),             (19)
```

again as large as `O_C(sqrt(R))`.  A bare geometry-of-numbers use of
the index cannot improve (12).

## 5. Audited conclusion

The primitive-chord transform contributes a reversible exact model
and the gcd-sensitive inequality (11).  It isolates a concrete
subproblem: any successful continuation along this route must exploit
the square condition (14) together with the congruences (5) in a way
that saves a power beyond their separate circle and lattice counts.
Convex strip width alone, the quadratic discriminant alone, and the
index of the transformed sublattice all retain an `R^(1/2)` loss when
the endpoint chord is primitive.

No radius-independent bound follows from this calculation.  In
particular it neither resolves the primitive-chord case nor proves the
global uniform `M(C)` conjecture.
