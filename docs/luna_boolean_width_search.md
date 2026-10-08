# Exact Boolean-width search (Luna audit)

This note records a bounded, integer-exact search for the Boolean width

\[
 F(A)=2^{-n}\sum_{x\in\{0,1\}^n}
       (\max_i a_i\!\cdot x-\min_i a_i\!\cdot x),
\]

where the rows are distinct integer vectors with coordinate sum one.  The
checker is [luna_boolean_width_search.py](/Users/cxc/Github/Jarnik/docs/luna_boolean_width_search.py).
It computes the integer numerator `W = 2**n F(A)` directly, and uses exact
integer elimination for affine dependence.  Because all rows have sum one,
linear dependence and affine dependence are equivalent: a linear relation
automatically has coefficient sum zero.

## Exhaustive box scans

For each `n` and bound `B`, rows are all points of
`[-B,B]^n ∩ {sum(a)=1}`.  The output below gives the minimum over all sets and
the minimum over dependent sets.  `W` is the sum of the widths over the
`2**n` masks; hence the displayed fraction is `F=W/2**n`.
The program also reports exact counts of violations of the general bound and
of the dependent-set bound; every count was zero in the scans below.

| n | B | rows | m | all minimum | dependent minimum |
|---:|---:|---:|---:|---|---|
| 2 | 3 | 6 | 2 | `W=2, F=1/2` | — |
| 2 | 3 | 6 | 3 | `W=4, F=1` | `W=4, F=1` |
| 2 | 3 | 6 | 4 | `W=6, F=3/2` | `W=6, F=3/2` |
| 3 | 4 | 60 | 2 | `W=4, F=1/2` | — |
| 3 | 4 | 60 | 3 | `W=6, F=3/4` | `W=8, F=1` |
| 3 | 4 | 60 | 4 | `W=8, F=1` | `W=8, F=1` |
| 4 | 1 | 16 | 2 | `W=8, F=1/2` | — |
| 4 | 1 | 16 | 3 | `W=12, F=3/4` | `W=16, F=1` |
| 4 | 1 | 16 | 4 | `W=14, F=7/8` | `W=16, F=1` |
| 4 | 1 | 16 | 5 | `W=18, F=9/8` | `W=18, F=9/8` |
| 4 | 2 | 80 | 2 | `W=8, F=1/2` | — |
| 4 | 2 | 80 | 3 | `W=12, F=3/4` | `W=16, F=1` |
| 4 | 2 | 80 | 4 | `W=14, F=7/8` | `W=16, F=1` |

The representative minimizing sets are retained in the program output.  For
example, at `(n,B,m)=(4,2,4)` the all-set minimum is attained by
`{(-2,-1,2,2), (-1,-2,2,2), (-1,-1,1,2), (-1,-1,2,1)}` and the dependent
minimum by `{(-2,-1,2,2), (-2,0,1,2), (-2,0,2,1), (-2,1,1,1)}`.

Commands used:

```text
python3 docs/luna_boolean_width_search.py --n 2 3 --bound 3 --m 2 3 4
python3 docs/luna_boolean_width_search.py --n 3 --bound 4 --m 2 3 4
python3 docs/luna_boolean_width_search.py --n 4 --bound 1 --m 2 3 4 5
python3 docs/luna_boolean_width_search.py --n 4 --bound 2 --m 2 3 4
```

No counterexample to either proposed inequality occurs in these exhaustive
scans.  For `n=3` with `m=2,3` and for `n=4` with `m=2,3,4`, the all-set
minima are exactly `1 - 2**(1-m)`; every dependent minimum is at least one.
The larger values in the `n=2`, `n=3,m=4`, and `m=5` rows reflect the limited
number of available independent directions or the extra row.

## Independent local audit of the distance-two case

The short-width proof reduces, after subtracting a common column baseline, to
equal-weight binary subsets.  The checker has an `--audit` mode that exhausts
all four-subsets of equal-weight binary rows in dimension four, retaining
those whose pair distance (positive mass of `a-b`) is at most two and whose
largest pair distance is two:

```text
python3 docs/luna_boolean_width_search.py --n 4 --bound 1 --m 4 --audit
```

There are two pair-distance multisets in this binary audit, and both have
minimum `W=16=2**4`.  Representatives of the two sharp patterns are

```text
{(1,1,0,0), (1,0,1,0), (1,0,0,1), (0,1,1,0)}
{(1,1,0,0), (1,0,1,0), (0,1,0,1), (0,0,1,1)}
```

Their six pair distances are respectively
`(1,1,1,1,1,2)` and `(1,1,1,1,2,2)`.  Direct mask enumeration gives `W=16`
and therefore `F=1` for each.  Exact ranks are 3 and 2, respectively, so
Direct exact affine-rank computation gives affine ranks 3 and 2,
respectively: the first is affinely
independent, while the second is dependent.  This independently checks the
two local four-row configurations used in the proof.  In particular, the
distance-two obstruction reaches width exactly one; it does not claim that
every sharp four-row pattern is dependent.  Adding either configuration to a
larger row set cannot reduce width.

For an additional proof audit, if all row-coordinate differences are at most
one, write each row as a fixed coordinatewise baseline plus the indicator of
an equal-size subset of coordinates.  If two subsets have Johnson distance
`s`, their pair width is exactly

```text
f(s) = s * binomial(2s,s) / 4**s,
```

so `f(1)=1/2`, `f(2)=3/4`, `f(3)=15/16`, and `f(4)=35/32`.  For three scalar
values, the range is half the sum of the three pairwise absolute differences.
Thus a distance-three pair plus any third distinct row has width at least
`(f(1)+f(2)+f(3))/2 = 35/32 > 1`; a distance-four pair already has width
`35/32`.  A distance-two pair forces every other row in a width-`<1` set to be
a common Johnson neighbor, and the audit above checks the two possible
four-row sharp patterns.  If all pair distances are one, the standard
Johnson-clique alternatives (a star or a costar) give the translated simplex
and width `1-2**(1-m)`.

The script uses only the Python standard library (`itertools`, `fractions`,
and `argparse`), so every value above is reproducible without floating-point
arithmetic or third-party packages.
