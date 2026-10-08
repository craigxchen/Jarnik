# Four point involutions and binary cut compatibility

This note records a local obstruction to obtaining clipped conductor from a
four point projective involution.  It uses the exact forward/inverse clipping
isometry in [the reciprocal stretch grid](mobius_reciprocal_stretch_grid.md).

Fix an odd split prime `p=pi bar(pi)`. Suppose the source and target
physical rows have binary `p`-allocations, with valuation widths `e` and
`e'`. After choosing actual anchors, write their signed phase valuations
as `t_i,t_i'`, with target extremes `L',U'`.
If a source valuation fiber contains two rows whose images
land at both target extremes, then the `p`-part of the clipped conductor is
trivial.

Indeed, the reciprocal clipping identity gives

```text
clip_I'(t_i') = F(clip_I(t_i)),
```

where `F` is an affine isometry and `I,I'` are the forward and inverse
coefficient intervals. The two source rows have the same `t_i`, so their
right sides agree.  Their target values are the two distinct extremes of
the target interval, hence `clip_I'(L')=clip_I'(U')`. A monotone interval
clipping map can do this only when the clipped target width is zero.  The
forward and inverse clipped widths are equal, so the forward width is zero
as well.

Conjugating the target if needed makes the determinant positive; this
reverses its two valuation classes without changing the conclusion.

Consequently, for a binary self-permutation with cut `S`, a necessary
condition for positive `p`-retention is that the permutation map each
source fiber wholly into one target fiber.  When source and target cuts are
the same, this says that the permutation preserves `S` or exchanges `S`
with its complement.

For a four-label double transposition this gives an immediate distinction.
With a `3|1` cut, every double transposition mixes the singleton with the
three-point fiber, so it forces zero clipped width at `p`. With a `2|2`
cut, a pairing internal to the two fibers preserves the cut and a pairing
with both pairs crossing complements it; these cases can retain conductor.
An all-same cut has no source `p`-width when these four rows form the
whole configuration. Four all-same seed rows inside a larger tuple do
not determine its clipped width.

## Exact fixture

Take `pi=2+i`, `p=5`, `e=2`, and auxiliary factors
`3+2i` and `4+i`. The pointwise primitive rows of squared radius `5525`

```text
z=(-14-73i, 62-41i, 22-71i, -14+73i)
```

have source half-angle rows, anchored at the first point,

```text
H=((1,0), (-3,-2), (-1,4), (73,14)).
```

Their signed `5`-valuations are `(0,0,0,2)`. The three exact trace-zero
matrices exchanging the prescribed pairs `(01)(23)`, `(02)(13)`, and
`(03)(12)` are

```text
(-6, 35; -4, 6),
(1, -38; -4, -1),
(-146, 193; -28, 146).
```

Their coefficient intervals at `5` are respectively `[-1,0]`, `[0,0]`,
and `[0,0]`.  Clipping the source interval `[0,2]` therefore gives width
zero in all three cases.  After reanchoring each target at its image of the
first row, the target signed valuations are respectively

```text
(0,0,2,0), (0,2,0,0), (0,-2,-2,-2).
```

Thus a `3|1` binary cut supplies a concrete local obstruction independent of
coefficient cancellation.  The companion `2|2` fixture has rows

```text
(-14-73i, 62-41i, 62+41i, -14+73i),
```

source valuations `(0,0,2,2)`, and the same pairings have clipped widths
`2,2,2`.  The exact calculations are in
[check_binary_four_involution_cut_compatibility.py](check_binary_four_involution_cut_compatibility.py).

This is a local compatibility result only.  It does not construct a global
map for an arbitrary endpoint tuple or prove the uniform point bound.
