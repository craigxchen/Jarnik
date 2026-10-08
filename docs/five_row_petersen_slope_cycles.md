# Cycle equations for the ten boundary slope ratios

These equations make the scale compatibility in the
[Petersen circulation dictionary](five_row_petersen_circulation.md)
explicit. The slope convention is fixed in the
[boundary restriction table](five_row_boundary_slope_relations.md).

Use the vertex order

```text
12, 13, 14, 15, 23, 24, 25, 34, 35, 45.
```

The following nine Petersen edges form a spanning tree:

```text
12-34, 12-35, 12-45,
13-24, 13-25, 13-45,
14-23, 14-25, 15-23.
```

The six chords are

```text
14-35, 15-24, 15-34, 23-45, 24-35, 25-34.
```

For a boundary vertex `S`, let its three incident edge coefficients be the
three entries of

```text
(alpha_S, beta_S, gamma_S),       gamma_S=-alpha_S-beta_S,
```

placed according to the outside-pair ordering: the edge whose disjoint pair
is the first two outside labels receives `gamma_S`, the edge for the first
and third receives `beta_S`, and the edge for the second and third receives
`alpha_S`.  Assume all three entries are nonzero and put
`r_S=alpha_S/beta_S`.

For an oriented cycle `C=(v_0,...,v_(m-1))`, define at each vertex

```text
rho_v = coefficient at the edge v-v_(next)
        / coefficient at the edge v-v_(previous).
```

If a nonzero vertex scale makes the two endpoint coefficients opposite on
every edge, telescoping gives the exact cycle equation

```text
product_(v in C) rho_v = (-1)^m.                       (1)
```

Conversely, fix a root scale and propagate vertex scales along the spanning
tree.  The equation for each chord is exactly (1) on its fundamental cycle.
If five chord equations hold, the remaining chord equation is automatic.
Write `s_v` for the propagated scale and `f_(v,e)` for the local edge
coefficient. The sum of all edge mismatches is

```text
sum_(e={v,w}) [s_v f_(v,e)+s_w f_(w,e)]
  = sum_v s_v sum_(e incident to v) f_(v,e) = 0.
```

Fourteen edge mismatches already vanish, so the last one vanishes too. Thus
five fundamental cycle equations are necessary and sufficient on the open
set where all three local coefficients are nonzero. The resulting
circulation determines a unique section up to common scale. For rational
slopes its coefficients are rational; clearing denominators gives an
integer section by the integral circulation identification. This is a
projective slope compatibility statement; it introduces no Gaussian arithmetic.

For the displayed tree, the six fundamental cycle lengths are `6,6,8,5,5,5`.
Omitting the third chord, `15-34`, gives five equations of degrees at most
`6,6,5,5,5` after clearing denominators. The sixth equation is
dependent, even though the Petersen cycle rank is six.

One explicit five-cycle is

```text
23 -> 14 -> 25 -> 13 -> 45 -> 23.
```

Its cleared slope equation is

```text
r_13 r_45 (r_23+1)(r_25+1) - r_23(r_14+1) = 0.         (2)
```

The [checker](check_five_row_petersen_slope_cycles.py) verifies (2) as an
exact polynomial identity in the six anticanonical `q`-coordinates,
rather than by numerical sampling.
