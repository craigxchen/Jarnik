# Ordered affine circuits retain the endpoint-versus-interior conductor

Let five distinct Gaussian integers `z_0,...,z_4` lie on one centered
circle, listed in their order on an arc. The primitive integer affine
relations on `0123` and `1234` span a sublattice of the full rank-two
affine-relation lattice. Write its index as `Q`. Every odd split-prime
allocation layer that separates the two endpoints `0,4` from all three
interior rows `1,2,3` divides `Q`, even after all accidental triangle
content is retained. After primitive Gaussian normalization, these layers
divide `gcd(Q,N)`, where `N` is the squared radius.

This is an exact evaluated-integrality statement. It applies to nested
prime-power allocations and arbitrary row units. It gives no upper bound
for `gcd(Q,N)` on an endpoint arc. Bounding the total index `Q` instead
would be too strong: an archived actual five-point endpoint family below
has `Q` unbounded but `gcd(Q,N)=1`.

The forced factor itself now has an
[exact subgroup-radius denominator formula](ordered_pair_triple_gap_denominator.md),
which removes accidental triangle valuations from the sufficient target.
The [accidental-conductor note](ordered_index_accidental_conductor.md)
classifies those additional valuations in the two-level case and gives
actual examples where the forced factor is smaller than `gcd(Q,N)`.

## Exact banded index

More generally, take `m>=4` distinct points on a circle. Write

```text
D_ijk=det(z_j-z_i,z_k-z_i),       g=gcd_(i<j<k)|D_ijk|,
T_j=D_(j,j+1,j+2),               0<=j<=m-3,
g_j=gcd_(I subset {j,j+1,j+2,j+3}, |I|=3)|D_I|,
                                      0<=j<=m-4.
```

All triangles are nonzero. Let `c^(j)` be the primitive circuit on the
four-row window beginning at `j`: its coefficient at `j+r` is
`(-1)^r D_(window minus {j+r})/g_j`, up to one common sign. The columns
`c^(0),...,c^(m-4)` form a basis over `Q` of the affine kernel. Their
integer span has exact index

```text
Q_band = g * product_(j=1)^(m-4) |T_j|
             / product_(j=0)^(m-4) g_j.                   (1)
```

Empty products are one. Indeed the first `m-3` rows of the circuit matrix
give a triangular minor whose absolute determinant is
`product_(j=0)^(m-4) |T_(j+1)|/g_j`. The
[complementary-minor identity](joint_affine_relation_lattice.md) makes
this minor equal to `Q_band |T_(m-3)|/g`, proving (1) and independence.
The opposite triangular minor also gives

```text
product_j |c^(j)_j|     = Q_band |T_(m-3)|/g,
product_j |c^(j)_(j+3)| = Q_band |T_0|/g.                 (2)
```

For five rows this specializes to

```text
Q = g |D_123|/(g_0123 g_1234).                           (3)
```

The local gcds must not be replaced by the global gcd. Equations (1)--(2)
identify exactly where products of local coefficient charges go: into
the nonsaturation index and a primitive triangle. In particular the
identity at any rational prime `p` is

```text
sum_j v_p(|c^(j)_j|) = v_p(Q_band)+v_p(|T_(m-3)|/g).      (4)
```

It does not supply an upper bound on that index. Dividing every point by
a common Gaussian factor divides every triangle and every displayed gcd
by the same integer norm, so all circuit vectors and indices are invariant.

## The nested conductor theorem

Fix an odd split prime `p=pi conjugate(pi)`, and put

```text
e=v_p(N),        t_i=v_pi(z_i),        0<=t_i<=e,
u=min(t_1,t_2,t_3),       v=max(t_1,t_2,t_3),
H_p=max(0,min(t_0,t_4)-v)+max(0,u-max(t_0,t_4)).          (5)
```

The integer `H_p` counts precisely the thresholds `h=1,...,e` for which
the cut `{i:t_i>=h}` restricts to `{0,4}` or `{1,2,3}`.
Then

```text
v_p(Q) >= H_p.                                          (6)
```

No two-level or squarefree hypothesis is imposed. If the Gaussian gcd
of the five points is removed, the new exponent of `p` in the squared
radius is `max_i t_i-min_i t_i`. All constant threshold layers disappear;
`H_p` is unchanged and is at most this exponent. Hence

```text
product_(odd split p) p^H_p divides gcd(Q,N_primitive).   (7)
```

A primitive Gaussian circle has no ramified or inert norm prime: such
a prime would divide every point. Nevertheless the local proof below
uses `p` odd and split; no assertion at two is silently included before
primitive normalization.

### Exact excess calculation

For a pair at equal valuation level, set
`r_ab=v_pi(z_a-z_b)-t_a>=0`. The exact triangle formula from
[the affine-radius note](affine_shape_radius_divisibility.md) is

```text
v_p(D_I)=e-range(t_i:i in I)
            +sum_(a<b in I, t_a=t_b) r_ab.               (8)
```

It follows from `2i D_ijk = +/- N (z_j-z_i)(z_k-z_i)
(z_j-z_k)/(z_i z_j z_k)`. At unequal levels the valuation of a
difference is their minimum; at equal levels its full excess is retained.
The unit `2i` is invertible at the chosen prime. Row units do not affect
the validity of this identity.

For a subset of at least three rows with at least three distinct
valuation levels, its triangle gcd has valuation exactly
`e-range`: choose a triple containing the minimum, maximum, and an
intermediate level. For a subset with exactly two levels, its triangle
gcd valuation is `e-range+min r_ab`, where the minimum is over its
equal-level pairs. A pair attaining the minimum, with a row at the
other level, attains this bound.

The range terms in (3) satisfy the elementary identity

```text
range(t_0,t_1,t_2,t_3)+range(t_1,t_2,t_3,t_4)
 -range(t_0,t_1,t_2,t_3,t_4)-range(t_1,t_2,t_3)=H_p.       (9)
```

If `H_p=0`, (6) follows from the integrality of `Q`. Suppose `H_p>0`.
Both endpoints are then strictly on one side of the entire interior
valuation interval.

If the interior occupies at least two levels, both four-row windows
and the full tuple occupy at least three levels. Their gcd excesses
are zero, and (3), (8), (9) give the exact formula

```text
v_p(Q)=H_p+sum_(1<=a<b<=3, t_a=t_b) r_ab >= H_p.          (10)
```

If the interior occupies one level, write its three pair excesses as
`r_12,r_13,r_23`, and let `a=min(r_12,r_13,r_23)`. Both window gcd
excesses equal `a`. The full-tuple excess is

```text
delta=0                         if t_0!=t_4,
delta=min(a,r_04)               if t_0=t_4.
```

Consequently the exact formula is

```text
v_p(Q)=H_p+r_12+r_13+r_23-2a+delta >= H_p.               (11)
```

This proves (6) including every accidental common factor. In the
two-level case with no excesses, the endpoint-versus-interior cut is
the only nonconstant five-row cut with a positive forced contribution.
Other cuts may still contribute through actual excesses.

## Singleton cuts charge actual relation coefficients

Suppose `sum_(i in J) c_i z_i=0`, with nonzero rational-integer
coefficients `c_i`. If one split prime has a singleton cut on `J`, isolating
row `i`, reduction at the Gaussian prime dividing all other rows forces
that Gaussian prime into `c_i`. Since `c_i` is a rational integer, the
full rational prime divides it. More generally a unique extreme valuation
gap `d` at row `i` forces `p^d | c_i`. Distinct prime charges multiply.
The rational-integer coefficient hypothesis matters: a Gaussian
coefficient would only receive the oriented Gaussian divisor.

For `b` repeated sign columns, let `n_i` count original labels whose
restriction to `J` isolates `i`, and put `h=|J|`. With one physical flip
per row, only the `h` flips at rows of `J` can change the restriction.
Thus at least `max(0,b n_i-h)` physical prime columns still isolate `i`.
If their prime log weights are at least `tau`, then

```text
log|c_i| >= tau max(0,b n_i-h).                          (12)
```

This improves the estimate that discards two entire copies of every
label. For local four-row circuits, multiplying their endpoint charges
without subtracting the exact index in (4) would overstate the resulting
bound on the saturated relation lattice.

## Two actual endpoint families distinguish the indices

### The rational quartic five-point family

The family in the [joint-relation note](joint_affine_relation_lattice.md)
has `R=Theta(n^8)` and every triangle polynomial of degree four. Its
global polynomial triangle gcd is constant. The gcd of triangle
polynomials on a four-subset has degree zero when row `0` is omitted,
and degree two when any other row is omitted. These claims are checked
by exact polynomial Euclidean division.

Therefore in the displayed label order `(0,1,2,3,4)`, the two local
gcds have degrees two and zero, and `Q=Theta(n^2)=Theta(R^(1/4))`.
That is not the geometric arc order. The distinct leading relative
phases in the [construction](four_row_rational_quartic_construction.md)
give the eventual arc order `(3,0,4,1,2)`. Both local gcds then have
degree two, and the middle triangle degree four, so `Q=Theta(1)`.

Here polynomial gcd degrees give integer-value growth rigorously:
divide by a primitive integral polynomial gcd and apply a fixed Bezout
identity to the remaining coprime polynomials. Their integer-value gcd
is bounded by a fixed nonzero integer. In particular the global gcd
being bounded does not justify assuming either local gcd is bounded.

The conductor part is not always one. Exact samples have
`Q=2937,gcd(Q,N)=89` at `n=12`, and `Q=261393,gcd(Q,N)=89` at `n=80`.
These are bounded-index fixtures, not a conjectured exact classification
of all parameter values.

### An ordered cubic-shape family has an unbounded external index

Use the [actual cubic-shape family](five_point_affine_shape_cubic_family.md)
with `T` a positive multiple of `600`. It is already Gaussian-primitive,
is ordered as `(0,1,2,3,4)` or its reversal on an arc shorter than
`5 sqrt(R)`, and satisfies `R~T^6/720`. Its triangle values are

```text
D_ijk=-(T/540) q_ijk(T),       g=T/15,
q_123=21T^2.
```

For the first window the four `q` polynomials are
`27T^2+108,56T^2+1008,50T^2+900,21T^2`; for the second they are
`21T^2,108T^2+108,112T^2+1008,25T^2+900`. Put `s=T/6`. Their exact gcds are

```text
a=108 if 3 does not divide s,       a=36 if 3 divides s,
b=36.
```

For verification, both gcds divide `756`, the gcd of the two-by-two
coefficient minors. Seven cannot divide the first gcd because it would
require `T^2=3 mod7`, or the second because it would require
`T^2=-1 mod7`. Their remaining two- and three-adic orders follow directly
after `T=6s` is substituted. Equation (3) now gives

```text
Q=7s^2  if 3 does not divide s,
Q=21s^2 if 3 divides s.                                 (13)
```

In particular `Q=Theta(R^(1/3))`. But write `T=60n`, with `10|n`.
The squared radius is

```text
N=product_(b in {60,30,20,15,12,10}) (1+b^2 n^2).
```

It is coprime to `n`, to two and five, and to the inert primes three
and seven. All prime factors of (13) lie among these primes. Thus

```text
gcd(Q,N)=1                                             (14)
```

for this entire actual endpoint family. Its normalized arc constant
tends to `2sqrt(5)`; it is not a counterexample for arbitrarily small
endpoint constants. More importantly, it does not refute a bound on
the conductor-supported index at any constant.

## What the positive conic metric charges

There is a direct circuit identity retaining the positive definite conic
metric. Use `P_1-P_0,P_2-P_0` as an affine coordinate basis. Divide its
integer Gram matrix by the gcd of its entries and write the result as
`[[A,B],[B,C]]`. The circle equation in these coordinates is

```text
A(u^2-u)+2Buv+C(v^2-v)=0.
```

For a primitive affine circuit `c_0,...,c_3` on the first four points,
the fourth point has coordinates `(-c_1/c_3,-c_2/c_3)`. Substitution and
`c_0+c_1+c_2+c_3=0` give the exact identity

```text
(A+C-2B)c_1 c_2 = -c_0(Ac_1+Cc_2).                    (15)
```

Thus `p^h|c_0` and `p` coprime to `c_1 c_2` force
`p^h|(A+C-2B)`, the primitive metric value of the chord `P_1P_2`.
Divisibility of the other circuit coefficients can weaken this inference.
It is not an additional factor of `p^h` in the metric content or its
determinant. For the first finite fixture below, the two circuits are

```text
(5,-17,13,-1,0),       (0,2,-3,6,-5),
[[A,B],[B,C]]=[[13,21],[21,34]],
AC-B^2=1,             A+C-2B=5=Q=gcd(Q,N).
```

The second fixture has circuits `(13,-16,8,-5,0)` and
`(0,4,-41,50,-13)`, primitive Gram matrix `[[1,9],[9,82]]`, determinant
one and chord value `65`. Its conductor prime thirteen also occurs in
that value exactly once. Both examples retain strict alternating circuit
signs and a positive definite primitive metric. They rule out charging
the conductor into the metric content or determinant on these hypotheses
alone; they do not exclude stronger quantitative use of the smaller arc
constant `1/2`. Equation (15) by itself gives no improvement to the
required conductor-index exponent.

## Finite arithmetic checks and remaining target

The [checker](check_ordered_affine_circuit_conductor_index.py) verifies
the general complementary-minor/index formula, singleton-prime charges,
and physical-count refinement on exact fixtures. It also checks the
nested local valuation formula on actual equal-norm Gaussian tuples,
including prime powers and accidental excess, and audits both endpoint
families above. A bounded scan of primitive short circle windows records
conductor overlap rather than treating ramified common scaling as new
overlap. Finite examples include

```text
N=1105, Q=5: (4,-33),(9,-32),(12,-31),(23,-24),(24,-23),
N=13325,Q=13: (-110,-35),(-109,-38),(-98,-61),(-94,-67),(-86,-77).
```

Their actual normalized arc constants are approximately `3.95544` and
`4.53598`; the checker certifies upper bounds four and `23/5` without
floating-point angle inequalities. These finite values prove neither
growth nor boundedness of the conductor part.

The unresolved target is an upper bound for `log gcd(Q,N)` on ordered
primitive five-point endpoint arcs in terms of their arc constant.
The divisibility theorem (7) supplies an actual prime-allocation lower
bound for that quantity. Multiplying local circuits without their
nonsaturation index, or bounding total `Q` instead, does not supply
the required upper bound. The separate
[ordered-sampling criterion](ordered_conductor_index_uniformity_criterion.md)
shows that even an estimate `log gcd(Q,N)<=delta log N+B` with fixed
`delta<1/15` on arcs of constant `1/2` would suffice for uniformity;
that arithmetic upper estimate remains unproved.
