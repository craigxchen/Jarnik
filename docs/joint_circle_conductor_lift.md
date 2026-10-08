# The joint involution on actual Gaussian circle directions

The six-point involution has exact degree-seven Gaussian representatives.
For its arc-preserving ordering their raw moduli are all `exp(64w+o(w))`.
Their actual least circle radius is at least `exp(30w-o(w))`, while the
old generic prime factors alone account for only `28w` in the logarithm.
The subsequent [standalone inflation theorem](joint_involution_standalone_inflation.md)
strengthens this to `exp((98/3)w-o(w))`, even when the original anchor is
removed and a common rotation is allowed. Every labeling of this
fixed-three pullback loses bounded normalized endpoint scale. The
identities and initial bounds leading to that theorem are recorded here.

## 1. Exact raw representatives

Use the primitive half-angle numerators `P_1,...,P_6`, their binary vectors,
and the real determinants `Delta_ij=det(v_i,v_j)`. The Segre forms `F,H,U,V,G`
are those of [joint_segre_involution.md](joint_segre_involution.md).
Normalize the first three rows to `infinity,0,1`. Pulling the transformed
coordinates `aH/V,bF/V,cH/F` back through this normalization and removing
the identically common factor `Delta_12` gives

```text
Q_4=H P_4-2 Delta_34 Delta_56 Delta_41 P_2,
Q_5=F P_5-2 Delta_35 Delta_46 Delta_51 P_2,
Q_6=H P_6+2 Delta_36 Delta_45 Delta_61 P_2.                (1)
```

Equivalently, for `(f_i,g_i)=(H,V),(F,V),(H,F)`,

```text
Delta_12 Q_i=Delta_i2 f_i P_1-Delta_i1 g_i P_2.           (2)
```

The identity `Delta_i2 P_1-Delta_i1 P_2=Delta_12 P_i`, together with
`A+B=Delta_12 Delta_35 Delta_46`, proves the cancellation in (1).
Each vector has total degree seven and degree two in its own moving row,
degree one in each other row. No numerical division remains in (1).

[check_joint_circle_lift.py](check_joint_circle_lift.py) verifies (1)--(2)
as identities in the twelve independent variables `P_i,bar(P_i)`. It uses
the bracket numerator `bar(P_i)P_j-P_i bar(P_j)`. The common denominator
`(2i)^3` is irrelevant at every odd core prime. Each Gaussian numerator
has only sixteen monomials.

## 2. Their size under the arc-preserving ordering

The [real-arc argument](joint_involution_real_arc_control.md) chooses one
of two opposite orders so that `F<0`,
and the six transformed directions lie in the original arc. In that order
`H,V` are positive sums of matching products, and have absolute size
`exp(48w+o(w))`. Formula (1) first gives

```text
|Q_i|<=exp(64w+o(w)),       |Im Q_i|<=exp(48w+o(w)).
```

The exact determinant identities from (2) give

```text
det(Q_i,P_1)=Delta_i1 g_i,
det(Q_i,P_2)=Delta_i2 f_i.                               (3)
```

For rows four and five use the first identity, whose form is `V`; for
row six use the second, whose form is `H`. Each right side has magnitude
`exp(64w+o(w))`. On the left, both directions lie in an angular interval
of width `exp(-16w+o(w))`, while `|P_j|=exp(16w+o(w))`. Consequently

```text
|Q_4|=|Q_5|=|Q_6|=exp(64w+o(w))                         (4)
```

in logarithmic scale. This lower bound uses the proved arc containment;
it does not presume absence of cancellation in the real part of (1).

## 3. The exact matching and Gaussian-content bridge

Put `W=(P_1,P_2,P_3,Q_4,Q_5,Q_6)`. Its five matching determinants equal
the quartic vector `T(M)` literally, not merely projectively. Indeed, in
the normalized chart the last three binary vectors are

```text
(aH,V),       (bF,V),       (cH,F).
```

Direct substitution gives `(-AFUG,-BUGV,CUGV,-DFGV,-EF^3)`.
Both sides have degree four in each original row and transform by
`det(A)^12` under a common invertible binary change `A`. Thus the chart
identity is the full polynomial identity.

Let `R_new` be the least Gaussian circle radius of these directions with
the original anchor retained, and let `z_0` be its integral anchor. Then
`z_i=z_0 W_i/bar(W_i)`. Each chord is

```text
z_i-z_j=2i z_0 det(W_j,W_i)/(bar(W_i)bar(W_j)),
```

with the harmless sign depending on the determinant convention. If
`g_T` is the ordinary integer gcd of the five entries of `T(M)`, the
Gaussian common factor of the fifteen three-chord matching products is

```text
gamma=unit*(2i)^3 z_0^3 g_T/product_i bar(W_i).           (5)
```

This is a nonzero Gaussian integer: the primitive integer matching vector
has an integer Bezout combination equal to one. Its logarithmic modulus
is exactly

```text
log|gamma|=3log R_new-sum_i log|W_i|+log g_T+3log2.       (6)
```

By (4), the sum in (6) is `240w+o(w)`. The actual local-content bound gives
`log g_T<=150w+o(w)`. Since `|gamma|>=1`,

```text
log R_new>=30w-o(w).                                    (7)
```

This requires no endpoint hypothesis on the output radius. It uses the
arc-preserving directions and their already proved raw sizes. Equivalently,
the output matching-height bound `42w-o(w)` and the actual angular width
`exp(-16w+o(w))` give the same inequality via three-chord products.

If these same directions were on an endpoint arc for their new radius,
the retained first three separations would also force
`log R_new<=32w+o(w)`. Then (6) would give

```text
log|gamma|<=6w+o(w).                                    (8)
```

Thus a universal lower bound exceeding `6w` for the retained old-prime
part of `gamma` would exclude such an endpoint output.

## 4. What the generic Gaussian divisor table actually says

For a core cut `S`, substitute a uniform prime order into each `P_i` with
`i in S`, while its formal conjugate has order zero. The minimal orders
of the expanded polynomials (1), and of their conjugates, give the generic
oriented valuations. The exact table is
[joint_circle_lift_cut_valuations.json](joint_circle_lift_cut_valuations.json).

For each moving row, summing its sixty-four generic orders gives

```text
sum ord_pi Q_i=64,       sum ord_barpi Q_i=32,
sum min(ord_pi Q_i,ord_barpi Q_i)=32.                    (9)
```

Thus its unavoidable generic rational content has logarithm `32w`, and
its total generic Gaussian content has logarithmic modulus `48w`.
The remaining Gaussian value can still have modulus `exp(16w+o(w))`.
Its prime factors and any further rational content have not been classified.

After generic conjugate-primitive reduction, every old output phase
allocation is zero or one. The total conductor width is fifty-six:
precisely the cuts meeting the unchanged rows `1,2,3` survive. The eight
cuts contained in `4,5,6`, including the empty cut, disappear. Hence the
generic old core contributes only `28w` to `log R_new`, exactly the least
conductor already required by the unchanged three rows and the anchor.

The corresponding six-point subset-size distribution gives generic old
common matching content `36w`. The sorted-valuation formula proving this
is [joint_circle_matching_content.md](joint_circle_matching_content.md).
Those generic phase valuations are not established as the actual ones:
extra cancellation in (1) can change either Gaussian orientation. In
particular, `36w` cannot yet replace the needed actual lower bound in (8).

## 5. Direct cuts settle the anchored and standalone cases

The initial `30w` lower radius bound left a narrow possible descent
window. The [five-cut inflation proof](joint_involution_five_cut_inflation.md)
closes that window without assuming the generic divisor table is actual.
The cuts `123,1234,1235,1236` each force at least `1.5w` of Gaussian
matching content, and `12` forces another `0.5w`, up to controlled errors.
This gives `log|gamma|>=6.5w-o(w)` and
`log R_new>=(193/6)w-o(w)` through (6).

Retaining the old anchor is essential to that particular five-cut
argument: common allocations of all six points can disappear when it is
deleted. The later [standalone theorem](joint_involution_standalone_inflation.md)
uses sixteen different cut bounds and exact complementary-cut symmetry
to prove `log|gamma|>=8w-o(w)` without that anchor. Consequently the
arc-preserving ordering requires `log R_new>=(98/3)w-o(w)`, or
`R_new>=R_old^(49/48-o(1))`. For every labeling of the fixed-three
pullback, its normalized arc length is at least `exp(w/3-o(w))`.

These bounds allow every common rotation and settle this specified
joint-map candidate. They do not address subsequent arbitrary
fractional-linear changes of its output directions, other joint maps,
or the general uniform circle-point bound.
