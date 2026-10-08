# Exact primitive content of the positive five-point Ptolemy system

Five cyclically ordered Gaussian-integer points give five coupled
positive rational cross-ratios. Their complete primitive cancellation
has the elementary description below. In particular, the three
endpoint-versus-interior matching numerators have one exact common
factor `F`, and the threshold gap norm `K` divides `F`.

This is a simultaneous content identity, not a growth theorem. The
direct radius bound for `F` has exponent `1/2`; the existing interior
triangle bound for `K` has the stronger exponent `1/4`. Neither gives
the sufficient exponent below `1/15` at arc constant `1/2`.

## 1. Positive cyclic conventions

Let `z_0,...,z_4` be distinct Gaussian integers on a circle centered at
zero, listed in cyclic order, and put

```text
Norm(z_i)=N=R^2,       l_ij=|z_i-z_j|,
u_i=l_(i,i+3) l_(i+1,i+2)/(l_(i,i+2) l_(i+1,i+3)),
```

with indices modulo five. Every `u_i` is in `(0,1)`. In each cyclic
four-point window, the denominator is the product of the diagonals;
Ptolemy says that it is the sum of the two positive matching products.
These ratios are rational even when the Gaussian row units differ:
the corresponding complex difference-product ratio is real, because
the two products have the same chord phase up to sign. It belongs to
`Q(i)`, hence to `Q`, and its absolute value is the displayed ratio.

Canceling common chord factors and then using Ptolemy gives

```text
u_(i-1) u_(i+1)=1-u_i.                              (1)
```

For example, `u_4 u_1=l_01 l_23/(l_02 l_13)=1-u_0`.
Thus the cyclic signs and the distinction between the diagonal and
the two boundary matchings are retained.

Write the first two ratios in lowest terms:

```text
u_0=a/b,       u_1=c/d,
0<a<b,        0<c<d,        gcd(a,b)=gcd(c,d)=1.
```

Equation (1) gives

```text
e=ad+bc-bd>0,
u_2=b(d-c)/(ad),
u_3=e/(ac),
u_4=d(b-a)/(bc).                                    (2)
```

Conversely these formulas satisfy all five equations (1) whenever
the stated inequalities and `e>0` hold. They describe the positive
rational pentagon system; they do not specify an integral circle
realization or its radius.

## 2. Exact simultaneous cancellation

Define

```text
H=gcd(b,d),       X=gcd(a,d-c),       Y=gcd(c,b-a).
```

Then `H,X,Y` are pairwise coprime and

```text
gcd(b(d-c),ad)=HX,
gcd(e,ac)=XY,
gcd(d(b-a),bc)=HY.                                  (3)
```

Moreover `HXY|e`, so put

```text
F=e/(HXY)>0.                                        (4)
```

Here is the complete list of primitive positive integer Ptolemy
triples. Each row means `numerator + complementary numerator =
denominator`; the three entries have gcd one.

| ratio | numerator | complementary numerator | denominator |
| --- | --- | --- | --- |
| `u_0` | `a` | `b-a` | `b` |
| `u_1` | `c` | `d-c` | `d` |
| `u_2` | `b(d-c)/(HX)` | `YF` | `ad/(HX)` |
| `u_3` | `HF` | `(b-a)(d-c)/(XY)` | `ac/(XY)` |
| `u_4` | `d(b-a)/(HY)` | `XF` | `bc/(HY)` |

For a positive rational number, let `num` mean its positive numerator
after reduction. The simultaneous content is therefore exactly

```text
F=gcd(num(1-u_2),num(u_3),num(1-u_4)).                (5)
```

### Primewise proof of (3)

Fix a rational prime `p` and write `v=v_p`. If `p|H`, then `p|b,d`
and `p` divides none of `a,c,d-c,b-a`. Thus `p` divides neither `X`
nor `Y`. If `p|X,Y`, then `p|a,c,d-c,b-a`, which would also imply
`p|b,d`, contradicting the two reduced fractions. This proves the
pairwise coprimality.

For the first gcd, if `p|b`, then `p` does not divide `a`. When
`p|d`, the number `d-c` is a unit, and its gcd valuation is
`min(v(b),v(d))=v(H)`. When `p` does not divide `d`, that valuation
is zero. If `p` does not divide `b` but divides `d`, the numerator
`b(d-c)` is a unit. In the remaining case both `b,d` are units, and
the gcd valuation is `min(v(a),v(d-c))=v(X)`. These cases prove the
first equality with full prime-power multiplicities. Exchanging
`(a,b)` with `(c,d)` proves the third equality.

For the middle gcd, a prime dividing both `a` and `c` does not divide
`e`, since `e=-bd mod p`. If `p|a` but `p` does not divide `c`, use
`e=ad-b(d-c)` and the fact that `b` is a unit. Truncating at `v(a)`
gives

```text
min(v(e),v(a))=min(v(d-c),v(a))=v(X).
```

This remains true when the two summands in `e` cancel to higher
order: the truncation is exactly what is needed. If `p|c` but not
`a`, the identity `e=bc-d(b-a)` similarly gives valuation `v(Y)`.
If neither `a` nor `c` is divisible by `p`, the gcd valuation is
zero. These exhaust the cases and prove the middle equality.

Finally, each of `H,X,Y` divides `e`, by its displayed linear
expressions. Their pairwise coprimality proves (4), and the table
and (5) follow.

## 3. The actual threshold gap divides the joint content

Now suppose the full Gaussian tuple is primitive. At an odd split
prime `p=pi bar(pi)` dividing `N`, put

```text
e_p=v_p(N),        t_i=v_pi(z_i),
U_O=[min(t_0,t_4),max(t_0,t_4)],
U_I=[min(t_1,t_2,t_3),max(t_1,t_2,t_3)],
h_p=distance(U_O,U_I),        K=product_p p^h_p.
```

This is the exact gap factor of
[the pair-triple denominator note](ordered_pair_triple_gap_denominator.md).
Full primitivity removes inert and ramified norm factors and gives
`min_i t_i=0,max_i t_i=e_p`.

For an equivalent formula requiring no chosen prime representatives,
let `A_O,A_I` be the norms of the Gaussian gcds of the endpoint pair
and interior triple, and put `N_O=N/A_O,N_I=N/A_I`. At each prime
their exponents are the interval lengths `r_O,r_I`. Since the convex
hull of the two intervals is `[0,e_p]`, their gap is
`max(0,e_p-r_O-r_I)`. Hence

```text
K=N/gcd(N,N_O N_I).
```

This is the representative-free formula used in the exact checker.

The three ratios in (5) are

```text
q_23=1-u_2=l_04 l_23/(l_03 l_24),
q_13=u_3  =l_04 l_13/(l_03 l_14),
q_12=1-u_4=l_04 l_12/(l_02 l_14).                    (6)
```

Each is the absolute value of the rational complex ratio with
numerator `(z_0-z_4)(z_a-z_b)` and denominator
`(z_0-z_b)(z_a-z_4)`, for its interior pair `a<b`.

Suppose `h_p>0`. If the endpoint valuation interval lies above the
interior interval, every denominator difference has unequal endpoint
valuations. Its valuation is exactly their minimum. The numerator
differences have valuations at least their corresponding minima.
Consequently

```text
v_p(q_ab)=v_pi(q_ab)
 >= min(t_0,t_4)-max(t_a,t_b) >= h_p.                (7)
```

If the endpoint interval lies below the interior interval, the same
calculation gives

```text
v_p(q_ab) >= min(t_a,t_b)-max(t_0,t_4) >= h_p.
```

Signs are units and do not affect these valuations. Any equal-level
excess is in a numerator difference and can only increase the result.
Thus `p^h_p` divides each of the three reduced numerators in (6).
Multiplying over the split primes and using (5) proves

```text
K | F.                                               (8)
```

This proof does not assume squarefree norm, common row units, or
absence of accidental chord content. It does not identify all prime
factors of `F` with the prescribed gap.

### Exact valuation at a gap prime

The excess over the prescribed gap also has a direct geometric
description.  Normalize the outer pair and the interior triple separately:

```text
G_O=gcd_G(z_0,z_4),             x_i=z_i/G_O       (i=0,4),
G_I=gcd_G(z_1,z_2,z_3),         y_i=z_i/G_I       (i=1,2,3).
```

Let `c_O` be the ordinary integer content of `x_4-x_0`, namely the gcd of
its two coordinates.  Let `c_I` be the ordinary integer gcd of the four
coordinates of `y_2-y_1` and `y_3-y_1`.  These definitions are independent
of the unit choices for the two Gaussian gcds.  At every source prime with
`h_p>0`, one has the exact identity

```text
v_p(F)=h_p+v_p(c_O)+v_p(c_I).                         (9)
```

Here the two accidental terms have a narrow scope.  The outer term can be
positive only when `t_0=t_4`, and the interior term can be positive only
when `t_1=t_2=t_3`.

To prove this, for an equal-level pair put

```text
rho_ij=v_pi(z_i-z_j)-t_i,
```

and put `rho_ij=0` when the two levels differ.  At unequal levels the
valuation of a difference is their minimum.  At equal levels the excess at
`bar(pi)` is the same `rho_ij`: after the common powers are removed, the two
local coordinates have equal product and are units, so inversion preserves
the valuation of their difference.

Suppose first that the outer interval lies above the interior interval.
Write

```text
A=min(t_0,t_4),       v=max(t_1,t_2,t_3),       A=v+h_p.
```

Every denominator pair in (6) crosses the gap, so it has no equal-level
excess.  For each of the three interior pairs `{a,b}` in
`{23,13,12}`, direct valuation of (6) gives

```text
v_p(q_ab)=h_p+[v-max(t_a,t_b)]+rho_04+rho_ab.          (10)
```

All these valuations are positive, so they are also the valuations of the
reduced positive numerators.  If the interior levels are not all equal,
one of the three pairs contains a row at level `v` and a row at a different
level.  For that pair both terms after `h_p+rho_04` vanish.  If
all three levels are equal, the minimum extra term is
`min(rho_12,rho_13,rho_23)`.  Reversing the allocation orientation proves
the same assertion when the outer interval lies below: in (10), replace
`v-max(t_a,t_b)` by `min(t_a,t_b)-min(t_1,t_2,t_3)`.

It remains to identify the two residue terms intrinsically.  After outer
normalization, unequal endpoint levels give normalized allocations zero
and `|t_4-t_0|`; their difference is not divisible by the rational prime
`p`.  At equal levels the rational content has `p`-valuation `rho_04`. Hence

```text
v_p(c_O)=rho_04 if t_0=t_4, and 0 otherwise.
```

For the normalized interior triple, a nonconstant allocation contains a
minimum-level row and a maximum-level row.  Their difference is not
divisible by rational `p`, and every pair difference is an integer linear
combination of the two differences defining `c_I`; hence `v_p(c_I)=0`.
For a constant allocation, the three pair differences have rational
contents of `p`-valuations `rho_ab`, and the same spanning observation gives

```text
v_p(c_I)=min(rho_12,rho_13,rho_23).
```

Taking the minimum of (10) proves (9), including arbitrary row units and
nested source powers.

The accidental part in (9) is not the accidental part of the ordered
affine circuit index `Q`.  The exact formulas in the
[affine-circuit note](ordered_affine_circuit_conductor_index.md) make the
difference explicit.  At a gap prime, if the interior levels are
nonconstant, then

```text
v_p(F)=h_p+v_p(c_O),
v_p(Q)=h_p+sum_(1<=a<b<=3, t_a=t_b) rho_ab.            (11)
```

If the interior levels are constant, their three ultrametric depths can be
written, after sorting, as `a,a,b` with `a<=b`.  Put
`r=v_p(c_O)`.  Then

```text
                         v_p(F)             v_p(Q)
t_0 != t_4              h_p+a              h_p+b
t_0  = t_4              h_p+a+r            h_p+b+min(a,r).   (12)
```

Thus `p^h_p` is forced into both invariants, but neither accidental excess
generally contains the other.  This comparison is local and exact; it is
not a new size estimate for either invariant.

## 4. Radius accounting

Suppose the ordered tuple lies on a minor arc of angular span `Delta`,
and let `Delta_I` be the span from row one to row three. Its arc length
is `R Delta`. Thus the endpoint hypothesis `R Delta<=C sqrt(R)` is
exactly `Delta<=C N^(-1/4)`.

Primitive norm is odd. All coordinates therefore satisfy `x+y=1 mod2`,
so every difference is divisible by `1+i`. The Gaussian numerator
and denominator products defining each `q_ab` are both divisible by
`(1+i)^2`. If `q_ab=n/d` is reduced, those products are, up to sign,
`n gamma,d gamma` for a Gaussian integer `gamma`. Bezout for `n,d`
shows that `gamma` is also divisible by `(1+i)^2`; hence `|gamma|>=2`.
It follows that

```text
num(q_ab) <= l_04 l_ab/2.
```

Now `l_04<=R Delta`, and one of the two successive interior chords
has length at most `R Delta_I/2`. Equation (5) therefore gives

```text
F <= N Delta Delta_I/4 <= (C^2/4) N^(1/2).            (13)
```

For comparison, let `A_I` be the norm of the Gaussian gcd of the
three interior points. Its prime exponent is
`v_p(A_I)=e_p-length(U_I)>=h_p`, so `K|A_I`. Dividing those
points by their gcd divides their triangle determinant by `A_I`;
the remaining nonzero determinant is an even integer, of absolute
value at least two. Meanwhile the chord-product formula and the two
successive angular gaps give

```text
|D_123|=l_12 l_13 l_23/(2R) <= N Delta_I^3/8.
```

Consequently the existing stronger estimate remains

```text
K <= A_I <= |D_123|/2
  <= N Delta_I^3/16 <= (C^3/16) N^(1/4).             (14)
```

At `C=1/2`, equations (13) and (14) read `F<=N^(1/2)/16` and
`K<=N^(1/4)/128`. The joint identity has not improved the exponent.
All norms here are the tuple's actual primitive norms. Dividing an
inherited tuple by a common Gaussian factor preserves its angles
and changes its normalized arc constant by `(N_new/N_old)^(1/4)`.

## 5. Actual circle fixtures and the separate shape family

The two actual primitive ordered tuples from
[the affine-circuit note](ordered_affine_circuit_conductor_index.md)
give:

| `N` | `u_0` | `u_1` | `(H,X,Y)` | `F` | `K` | affine index `Q` |
| --- | --- | --- | --- | --- | --- | --- |
| `1105` | `1/2` | `51/52` | `(2,1,1)` | `25` | `5` | `5` |
| `13325` | `40/41` | `1/2` | `(1,1,1)` | `39` | `13` | `13` |

The tuples are respectively

```text
(4,-33),(9,-32),(12,-31),(23,-24),(24,-23),
(-110,-35),(-109,-38),(-98,-61),(-94,-67),(-86,-77).
```

The first has `|D_123|=10`, so `F=25` cannot be substituted for `K`
in the triangle estimate `K<=|D_123|/2`. Their normalized arc constants
are approximately `3.95544` and `4.53598`; they are outside `C=1/2`.

The first fixture also realizes the distinction in (11).  At
`pi=2+i` above five its allocation is

```text
(t_0,...,t_4)=(1,0,0,0,1),       h_5=1.
```

After separate normalization its outer pair is
`(-5-14i,5-14i)`, so `c_O=10`, while `c_I=1`.  Therefore (9) gives
`v_5(F)=1+1=2`.  The affine index has only `v_5(Q)=1`: the extra outer
chord content belongs to `F`, not to `Q`.

Separately, the positive rational pentagons

```text
a=c=3m+1,       b=d=4m+1,       m>=1
```

have `H=4m+1,X=Y=1,F=2m+1`, while

```text
(u_0,u_1,u_2,u_3,u_4) -> (3/4,3/4,1/3,8/9,1/3).
```

Thus joint primitive content can grow while every ratio stays away
from zero and one. This is a shape-only family. No integral circle
realization, radius estimate, or endpoint arc assertion is made for it.

## 6. Relation to the existing routes and verification

[Ordered primitive chord residues](ordered_residue_growth.md) retain
each four-row positive Ptolemy identity; the present calculation
normalizes all five simultaneously and identifies their exact shared
numerator content. The
[nonlinear consecutive audit](ordered_consecutive_ptolemy_nonlinear_audit.md)
shows why a bounded shape surplus does not itself give radius growth.
The [cross-ratio height norm](coupled_crossratio_height_norm.md) already
retains the universal projective exchanges and their cut valuations.
Formula (5) is exact evaluated gcd bookkeeping within that same system,
not an additional independent equation or an exponent gain.

The missing statement is still a radius-sensitive estimate on the
gap-supported part: `log K<=delta log N+B` with `delta<1/15` on
primitive ordered arcs of constant `1/2`. Positivity and (3)--(8)
alone have not supplied it.

The accompanying
[exact checker](check_positive_pentagon_primitive_content.py) verifies
the cancellations over 605,550 positive reduced rational pairs, all
five primitive triples, the exact circle fixtures and their Gaussian
gcd gap factors, the exact excess formula (9) on four literal tuples with
varying units and nested source powers, the parity-based numerator bounds,
and the separate shape family. The primewise proofs above establish the
general claims; finite checks supplement them.
