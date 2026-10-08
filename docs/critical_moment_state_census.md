# Exact size of the compressed moment state space

For each fixed parity and row-group split in the
[compressed automaton](critical_moment_compressed_automaton.md), there
are at most **326,496 valid states**, substantially fewer than the coarse
bound `1215*8!`. This counts every valid state, whether reachable or not.
It does not prove a weighted-path exclusion or a lattice-circle bound.

The exact counts are the same for both parities of `s`:

| Number `q` of positive cross orientations | Valid height vectors | Valid states |
|---|---:|---:|
| 0 or 4 | 255 | 322,560 |
| 1 or 3 | 262 | 325,584 |
| 2 | 264 | 326,496 |

Here is a short exhaustive counting proof. There are only fifteen possible
height vectors on the four `A` rows, because `r_0=0` and their diameter
is at most one. Each of the four `B` heights belongs to
`{0,delta_b,2delta_b}`. Enumerate those `15*81=1215` possibilities.

For each pair of rows, the two possible comparison signs specify the two
possible disagreement bits `e=0,1`. Reconstruct the pair phase by

```text
k_ij = D_ij+r_i+r_j-delta_i-delta_j+2e  (mod 4).
```

Reject a comparison unless it satisfies the corresponding critical or
successor height table and places negative normalized values before
positive ones. If neither comparison survives, the height vector gives
no state. If only one survives, impose that directed precedence relation.
If both survive, impose no relation. Thus the compatible permutations
are exactly the linear extensions of these precedence constraints.
This equivalence uses every local condition in the original validator;
it makes no assumption about transitivity of the unforced comparisons.

For a subset `U` of the eight rows let `f(U)` count permissible orderings
of exactly `U`. Set `f(empty)=1`. For each row `i` outside `U` whose
predecessors all lie in `U`, add `f(U)` to `f(U union {i})`. Induction on
`|U|` proves that `f(all rows)` counts every compatible permutation
exactly once. Cyclic constraints automatically contribute zero. The
entire calculation uses only 256 subsets for each height vector.

The histogram of nonzero extension counts for `q=0,4` is

```text
extensions per height:  576  720  1440  5040  40320
number of heights:       70  112    56    16      1
```

For `q=1,3`, add one height with 144 extensions, three with 240 and
three with 720. For `q=2`, add one with 96, four with 240 and four
with 720. Summing these exact histograms gives the displayed table.

The [checker](check_critical_moment_state_census.py) reproduces all ten
cases by integer subset dynamic programming. It additionally compares
the predecessor test with the original state validator on 97,200
deterministically sampled height/permutation pairs. These samples audit
the implementation; the complete counts come from the exhaustive
height enumeration and the proven linear-extension recurrence.
