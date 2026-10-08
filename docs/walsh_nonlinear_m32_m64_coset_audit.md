# Finite nonlinear Walsh coset gain audit at 32 and 64 rows

This note records a bounded, fixed-seed computation for the affine-coset
character statistic

```text
gain(P,l) = 2 k(P,l) - h,
 k(P,l) = #{x in P : a(x)|V = l},   h = |V|.
```

For a repeated physical profile with `b` copies of each nonzero Walsh label,
the exact physical coefficient sum is

```text
L(P,l) = b M/2 + h - 2 k(P,l).
```

Thus `gain` is the reduction relative to the unaligned value `b M/2`.
The computation below concerns this statistic only; it makes no claim about
endpoint arcs or a general nonlinear assignment theorem.

## The supplied M=32 map

The map is the assignment in
[the preceding nonlinear obstruction note](walsh_nonlinear_coset_height_obstruction.md),
using every label `1,...,31` once and label `3` twice. All subspaces and all
cosets are enumerated. For each coset, the best gain is the maximum over all
nonzero restrictions `l`.

| `h` | directions | cosets | maximum `k` | maximum gain | mean best gain |
|---:|---:|---:|---:|---:|---:|
| 4 | 155 | 1,240 | 3 | 2 | `-213/620` |
| 8 | 155 | 620 | 4 | 0 | `-111/31` |
| 16 | 31 | 62 | 3 | -10 | `-368/31` |

The complete best-gain histograms are

```text
h=4:  -2:341,  0:771,  2:128
h=8:  -6:17,  -4:458, -2:143, 0:2
h=16: -12:58, -10:4
```

In particular, the 128 positive `h=4` cosets are exactly the 128 plane
triples with `k=3` in the original checker. There are no positive cosets at
`h=8` or `h=16`.

## Fixed-seed M=64 examples

Each assignment uses `[1,1,2,...,63]`, so it is balanced in the stated sense.
There are 651 dimension-2 directions, with 16 cosets per direction, and
1,395 dimension-3 directions, with 8 cosets per direction. The total numbers
of cosets are therefore 10,416 and 11,160.

For a uniformly shuffled assignment with seed `20260915`, the results are

| dimension | positive cosets | positive fraction | best-gain histogram |
|---:|---:|---:|---|
| 2 (`h=4`) | 1,530 | 14.69% | `-4:30, -2:2849, 0:6007, 2:1422, 4:108` |
| 3 (`h=8`) | 50 | 0.448% | `-6:215, -4:6607, -2:3779, 0:509, 2:50` |

Independent random seeds `20260916` and `20260917` gave positive counts
`(1452,31)` and `(1568,50)` for dimensions `(2,3)`.

For an affine comparison, seed `20260918` generated the full-rank columns
`[57,20,56,10,40,59]` and offset `55`; the unique zero output was replaced by
the duplicate label `1`. Disjoint row swaps preserve the exact balanced label
multiset. Using swap seed `2026091800 + s` for `s` swaps gives

| disjoint swaps | positive `h=4` cosets | positive `h=8` cosets |
|---:|---:|---:|
| 0 | 312 | 10 |
| 1 | 446 | 15 |
| 2 | 582 | 28 |
| 4 | 728 | 21 |
| 8 | 1,049 | 28 |
| 16 | 1,439 | 42 |

The corresponding exact best-gain means are

```text
swaps  h=4             h=8
  0    -967/1302       -6266/1395
  1    -3551/5208      -4775/1116
  2    -157/248        -3823/930
  4    -69/124         -21613/5580
  8    -191/434        -19741/5580
 16    -359/1302       -4481/1395
```

The dimension-2 positive counts grow to roughly 14% after 16 swaps, while
dimension-3 positives remain below 0.4% in these examples. These are sample
statistics, not a density statement for all balanced assignments.

## Restriction-count balance

The positive `h=4` coset count vectors, sorted and with one entry for each of
the three nonzero restrictions, are only

```text
(3,0,0), (3,1,0), (4,0,0).
```

For `h=8`, positive vectors have one entry at least 5; examples from the
enumerations include `(5,1,0,0,0,0,0)`, `(5,2,1,0,0,0,0)`, and
`(6,0,0,0,0,0,0)`. Hence no positive coset in these computations has exact
uniform nonzero-restriction counts. The simple divisibility check already
rules out exact per-coset equality (`4` is not divisible by `3`, and `8` is
not divisible by `7`), while the enumerated vectors show substantial rather
than near-uniform concentration.

As a separate direct tally, each positive coset has a unique winning
restriction in dimensions 2 and 3. Counting those winners direction by
direction, exact equality across all restrictions occurred for 60/651
dimension-2 directions in the random seed-20260915 assignment, 17/651 for
the unperturbed affine assignment, and 73/651 after 16 swaps. It occurred for
0/1,395 dimension-3 directions in all three cases. This is only an exact
finite tally of the unique winners; it is not an LP feasibility result.

## Reproduction

Run

```text
python3 docs/check_walsh_nonlinear_m32_m64_coset_audit.py
```

The checker uses exact integer counts and `fractions.Fraction` for every
reported mean. It enumerates dimensions 2--4 at `M=32`, and dimensions 2--3
at `M=64`; the largest search has `1,395 * 8 * 7 = 78,120`
coset/restriction triples per assignment.
