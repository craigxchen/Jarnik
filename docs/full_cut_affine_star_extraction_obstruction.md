# Full cuts prevent affine extraction of a nontrivial star symmetry

Fix two output halves `A,B`, each of size `k>=2`. Consider integer
multiplicative affine combinations of a source tuple whose every
nonconstant cut class has a nonunit Gaussian block, with disjoint rational
split-prime supports between different classes. Every coefficient row has
sum one. If the normalized outputs have an exact half-symmetric threshold
profile, then either the coefficient rows repeat within each half or the
output profile has no star layer at all.

Consequently this operation cannot produce `2k` distinct Gaussian integer
outputs on an arc of length `C sqrt(R)` with arbitrary row units and
`C<=sqrt(2)`, or with one literal common row unit and `C<=2`.
The statement holds for every `k`, not just four outputs. It restricts a
proposed extraction operation; it proves no general lattice-point bound.

The hypothesis concerns the complete source blocks. A fair leading core
with uncontrolled correction factors, or different source cut layers at
the same rational prime, is not covered by this argument.

For a literal disjoint full fair core with small Gaussian row factors,
[the correction-height bound](full_fair_half_product_correction_height.md)
shows that the reduced half-product correction left after removing the
star has positive radius-scale logarithmic norm. Thus the core extraction
does not supply the `o(W)` correction needed by a small-height
perturbation of exact half symmetry.

## 1. The exact condition on a valuation column

For a real vector `v` indexed by `A union B`, call it regular when the two
multisets of values agree:

```text
{v_i : i in A} = {v_i : i in B}, with multiplicity.       (1)
```

Call it a star vector when it is constant on each half. Constants belong
to both classes; a nonconstant star has distinct values `alpha,beta` on
the two halves. Let `H` be the union of these two classes.

All nonconstant threshold cuts of `v` have either equal counts in the two
halves or are precisely `A` or `B` if and only if `v in H`. Indeed, if all
thresholds have equal counts, the finite value distributions agree. If a
star threshold occurs, every other nontrivial threshold is comparable to
it. A proper nonempty subset of one half, or a proper superset of that
half, cannot have equal counts in the two halves. Thus there can be no
other distinct nontrivial threshold: `v` is constant on each half.

This criterion is unchanged by adding a constant vector, changing sign,
or multiplying by a positive scalar. It therefore survives the exact
primewise common-gcd normalization.

## 2. A nonconstant star forces every coefficient column to be a star

Let `a_1,...,a_m` be the output vectors of coefficients of the individual
source rows. For a source subset `T`, write

```text
n(T)=sum_(j in T) a_j,             n(all)=1.              (2)
```

Disjoint source prime supports give residual valuations proportional to
`n(T)-min n(T)` and `max n(T)-n(T)` at the block for `T`. Thus exact half
symmetry implies `n(T) in H` for every subset. Only one orientation of
each source cut is needed, since `n(T^c)=1-n(T)`. The empty and full
subsets give constants automatically.

Suppose some `n(T)=s` is a nonconstant star with values `alpha!=beta`.
For each source label `j`, both `a_j` and

```text
s+a_j=n(T union {j})      if j is outside T,
s-a_j=n(T minus {j})      if j is inside T                (3)
```

belong to `H`. We claim that `a_j` must be constant on each half.
Otherwise it is regular, with the same non-singleton finite multiset `U`
on the two halves. The two half-value multisets in (3) are
`alpha+U,beta+U`, or `alpha-U,beta-U`. They cannot agree: a nonzero
translation does not preserve a finite nonempty multiset (compare their
minima). Nor is (3) a star, since `U` has at least two distinct values.
This contradicts membership in `H`.

Every `a_j` is therefore constant within each half. The coefficient rows
within either half are identical, so their formal Gaussian products are
equal, including any source units. Clearing a common denominator and
dividing a common Gaussian gcd preserves those repetitions. This proves
the asserted alternative.

The same combinatorial proof works for any common row degree, but degree
one is the radius-preserving affine operation considered here.

## 3. A regular-only output cannot satisfy the endpoint pair inequalities

Suppose the outputs are distinct, so Section 2 excludes every star layer.
Retain their complete common Gaussian divisor `d`. Set

```text
W=log R^2,        D=log Norm(d),        W_var=W-D.
```

A regular threshold has `t` selected rows in each half. Among the `k^2`
cross-half pairs it separates exactly `2t(k-t)`, at most `k^2/2`.
Summing with the actual prime-layer weights gives

```text
(1/k^2) sum_(i in A,j in B) log P_ij <= W_var/2 <= W/2.  (4)
```

But the exact chord bound in
[linear_allocation_affine_rigidity.md](linear_allocation_affine_rigidity.md)
gives `log P_ij>W/2` for every pair on the stated short arcs. Strictness
also holds at the boundary values of `C` because a nonzero chord is
strictly shorter than its containing arc. Equation (4) is a contradiction.
No phase approximation theorem is needed for this last step.

## 4. Why the no-star alternative is necessary for more than four outputs

There are six distinct coefficient rows of sum one for which every
source subset sum is regular. With four source columns, take

```text
A: (-2,1,2,0), (-1,2,0,0), (0,0,1,0),
B: (-1,0,2,0), (0,1,0,0), (-2,2,1,0).                  (5)
```

The first three coordinate multisets agree on the two halves. Each row
has sum one, so every two-coordinate sum is one minus the remaining
coordinate. The fourth column is zero. Thus all sixteen subset sums are
regular. With disjoint nonunit full-cut source blocks the normalized
products are distinct, since the singleton source classes detect every
distinct coefficient row.

This example shows why the four-output obstruction in
[affine_full_fair_to_half_symmetric_obstruction.md](affine_full_fair_to_half_symmetric_obstruction.md)
does not extend to a blanket impossibility of all half-symmetric products.
What is impossible in every size is a nonconstant star with distinct
outputs; regular-only products fail the endpoint inequalities instead.

## Verification

The standard-library checker
[check_full_cut_affine_star_extraction_obstruction.py](check_full_cut_affine_star_extraction_obstruction.py)
checks the threshold characterization and the shifted-multiset argument,
then realizes (5) using seven distinct rational split primes. It verifies
the primitive equal-norm Gaussian outputs and their cross-pair norm
budget in exact integer arithmetic. These checks verify the identities
and the scoped example, not a general endpoint theorem.
