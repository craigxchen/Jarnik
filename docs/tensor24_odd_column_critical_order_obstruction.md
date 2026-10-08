# The fixed tensor code fails the critical odd-moment sign order

The eight-row, 22-column sign matrix in
`tensor24_odd_column_certificate_gap_fixture.json` has no nonzero real
coefficients `a_j` with pairwise distinct absolute values satisfying

```text
S(a^(2k+1)) = alpha_k 1,       0 <= k <= 4.               (1)
```

The exclusion uses the [critical odd-moment root-order lemma](critical_odd_moment_root_order.md),
which strengthens the sign-variation test already passed by this code.
It concerns this fixed sign template. It makes no claim about other codes,
endpoint arcs, or a uniform lattice-circle bound.

## Sixteen simultaneous necessary patterns

For each pair of rows `i,j` differing in exactly eleven retained columns,
put `c=(S_i-S_j)/2`. Its eleven nonzero entries are `+1` or `-1`, and (1)
gives

```text
sum_(k:c_k!=0) (c_k a_k)^(2r+1) = 0,    0 <= r <= 4.
```

The equality holds because odd powers preserve the sign `c_k`. Thus, when
the eleven supported columns are read in increasing order of `|a_k|`, the
signs `c_k sign(a_k)` must be either

```text
P = ++--++--++-
```

or `-P`. The matrix has sixteen such row pairs. They are exactly the
cross-pairs between row sets `{0,2,3,4}` and `{1,5,6,7}`. Every pair
must use the same global order of all 22 magnitudes and the same 22
coefficient signs.

## Exhaustive finite obstruction

The [exact checker](check_tensor24_odd_column_critical_order.py) searches
every global column order and coefficient-sign assignment through prefix
states. A state records the selected-column mask and, for each row-pair
certificate already touched, whether its allowed pattern is `P` or `-P`.
The position within that pattern is the number of previously selected
columns in its support. At each step the checker tries every unselected
column and every sign consistent with all touched certificates; a first
occurrence fixes that certificate's orientation.

This state description is complete: two prefixes with the same mask and
orientation bits have identical allowed continuations. Memoization is
therefore exact. The search exhausts 3,029 states and reaches no full
22-column assignment. No arithmetic approximation enters the check.

Consequently (1) is impossible under nonzero, pairwise distinct absolute
coefficient values. This resolves feasibility for this fixed code without
using the distance-twelve pairs or any numerical solve.

A [separate ternary-state implementation](tensor24_odd_column_critical_order_core.md)
reproduces the full search and proves that ten of the sixteen certificates
already suffice. The [family extension](tensor24_odd_column_family_critical_order.md)
excludes all 21,120 balanced deletions in the normalized odd-positive tensor
family by reducing them to four equivalence classes.
