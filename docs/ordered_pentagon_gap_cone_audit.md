# Ordered pentagon gap inequalities miss the balanced cut profile

## Status

This is a restricted inequality-cone audit.  Its generators are the proved
five-row gap estimates

```text
log K_P <= W/4-log A_P+O_C(1),
log K_P <= log F_P <= W/2-log A_P+O_C(1),            (1)
```

for ordered five-subsets `P`, where `A_P` is the norm of the common Gaussian
divisor of the retained rows and `W_P=W-log A_P` is their primitive squared-
radius logarithm.  These formulas retain the improvement in the normalized
arc constant after division.  Arbitrary nonnegative
linear combinations are allowed, including selections whose interior
triangles are disjoint.  On the complete equally weighted cut profile, every
triangle generator uses more than three times as much radius exponent as the
gap mass it detects; the `F` generator is worse.  Moreover many cut classes,
including every singleton cut, are invisible to the whole cone.

Thus multiplication or weighted selection of the existing pentagon content
and interior-triangle bounds cannot yield a cardinality penalty on the
balanced profile.  This does not rule out inequalities using signed
cancellation, new coupled residues, or a quantity that is nonzero on the
missing cuts.

## 1. Cut incidence of one ordered pentagon

Fix `m>=5` linearly ordered rows.  Represent an unoriented nonconstant cut
by the side `S` not containing row zero.  There are

```text
D_m=2^(m-1)-1                                             (2)
```

such cuts.  Give every cut a core block of logarithmic norm `w`, so
`W=D_m w`. This is the idealized leading valuation profile in the uniform
profile extraction. No exact equal-block Gaussian realization or endpoint
localization is asserted here. For fixed `m`, replacing the block weights
by `w+o(w)` changes each of the finitely many quantities below by `o(w)`;
the strict exponent deficit therefore persists in that extracted regime.

For an ordered five-subset

```text
P={p_0<p_1<p_2<p_3<p_4},
O_P={p_0,p_4},             I_P={p_1,p_2,p_3},
```

its gap has full exponent `w` at a binary cut precisely when `O_P` lies on
one side and `I_P` on the other.  The other `m-5` rows are free.  After
identifying a cut with its complement, the number detected is exactly

```text
2^(m-5).                                                  (3)
```

Consequently, on the balanced profile,

```text
log K_P=2^(m-5)w=[2^(m-5)/D_m]W.                         (4)
```

The primitive radius of the retained five rows discards precisely the cuts
that are constant on those rows.  Each of the `2^4-1=15` nonconstant
unoriented cuts of the five labels has `2^(m-5)` extensions.  Hence

```text
W_P=15*2^(m-5)w.                                        (5)
```

Put `s=2^(m-5)`.  Since `D_m=16s-1`, the discarded constant-cut mass is

```text
log A_P=W-W_P=(s-1)w.                                   (6)
```

Common-divisor removal changes the normalized arc constant to
`C_P=C A_P^(-1/4)`.  Substitution in the actual estimates
`K_P<=C_P^3 N_P^(1/4)/16` and
`F_P<=C_P^2 N_P^(1/2)/4` gives (1).  On the balanced profile their leading
costs are therefore

```text
triangle: (3s+3/4)w,             F: (7s+1/2)w.           (7)
```

This is an exact valuation count.  Equal-level chord excesses may enlarge
`F_P`, but do not enlarge the prescribed gap `K_P` and hence cannot improve
the left side of (1).

## 2. Every nonnegative combination has the same deficit

Give the triangle inequalities in (1) arbitrary weights `lambda_P>=0`.
Equations (1) and (4) evaluate on the balanced profile as

```text
sum_P lambda_P log K_P
  = s w sum_P lambda_P,

sum_P lambda_P (W/4-log A_P+O_C(1))
  = (3s+3/4)w sum_P lambda_P+O_C(sum_P lambda_P).        (8)
```

The leading ratio of detected gap mass to the available upper bound is

```text
alpha_m=4s/(12s+3)<1/3.                                  (9)
```

The ratio is independent of the weights and of overlap among the interior
triangles.  Restricting to disjoint triangles merely deletes generators and
does not change (9).  For `m=6,7` it is respectively `8/27` and `16/51`.
The `F` efficiency is `2s/(14s+1)`, which is smaller.

Equivalently, even if all nonconstant cuts were coverable, an averaged
fractional cover would need total generator weight at least
`D_m/2^(m-5)`, and its triangle right side would be at least

```text
[(3s+3/4)/s]D_m w=(3+3/(4s))W,
```

before any loss from ordered cuts that cannot be covered at all.

## 3. The actual ordered cone has blind cuts

The preceding average is generous: the ordered generators do not cover all
cuts.  Every detected cut contains at least the two rows of `O_P` on one
side and the three rows of `I_P` on the other.  Hence every singleton cut is
invisible for every `m`.  Endpoint ordering creates further blind pair and
larger cuts.  In fact a binary cut word is blind exactly when it contains
neither `01110` nor `10001` as a subsequence.  Such a subsequence is exactly
a five-subset whose ordered endpoints have one bit and whose three interior
rows have the other.

There are exactly

```text
10m-36                                                    (10)
```

blind unoriented nonconstant cuts for every `m>=5`.  Count oriented words
first.  If both endpoint bits are zero, a blind nonconstant word has one or
two ones.  With `n=m-2` interior positions, there are `n` choices for one
one.  Two ones must be at distance at most three, or they form `10001`; the
count is

```text
(n-1)+(n-2)+(n-3)=3n-6.
```

Thus there are `4m-14` blind words with zero endpoints, and the same number
with one endpoints.

For opposite endpoints, take the zero-to-one orientation and remove the
maximal initial zero-run and final one-run.  The remaining core is empty, or
starts in one, ends in zero, and contains at most two of each bit.  Its only
possibilities are

```text
empty, 10, 100, 110, 1100, 1010.
```

Distributing the two nonempty endpoint runs gives

```text
(m-1)+(m-3)+2(m-4)+2(m-5)=6m-22
```

words; the last term is zero at `m=5`.  The reverse orientation contributes
the same number.  Altogether there are `20m-72` oriented nonconstant blind
words.  Complementation acts freely, proving (10).

Direct enumeration gives:

| rows | ordered five-subsets | cut classes | incidences per generator | blind cuts |
| ---: | ---: | ---: | ---: | ---: |
| 6 | 6 | 31 | 2 | 24 |
| 7 | 21 | 63 | 4 | 34 |

Thus no nonnegative combination of these gap functionals is coefficientwise
positive on every nonnegative cut weight: it vanishes on the blind rays of
the cut cone.  On the fixed balanced profile, (8)--(9) separately show the
strict radius-exponent deficit.  Together these facts prevent a contradiction
from multiplying or selecting the listed five-window inequalities.

The blind proportion `(10m-36)/(2^(m-1)-1)` tends to zero exponentially.
Blind support is therefore a finite-row coefficientwise obstruction, not a
positive-density obstruction for large `m`.  The radius-cost deficit in
(8)--(9) is separate and persists as `m` grows.

The enumeration is checked by
`check_ordered_pentagon_gap_cone.py`.  It uses only cut incidence, not a
probabilistic approximation.

## 4. Scope of the obstruction

This calculation applies to inequalities whose left sides are nonnegative
combinations of the prescribed five-row endpoint-versus-interior gaps and
whose right sides are exactly the two inherited-geometry estimates in (1).
Pair, Vandermonde, or other inequalities that lower these radius costs are
outside the audited cone and require a separate analysis; no claim about
such an enlarged cone is made here.

The adjacent-pentagon positive exchange also does not add such a generator:
on common gap support its primitive common numerator cancels both prescribed
gap powers.  A successful global budget therefore needs a quantity with
different cut support, a signed relation that exploits cancellation rather
than discarding it, or a genuinely smaller radius exponent per detected cut.
