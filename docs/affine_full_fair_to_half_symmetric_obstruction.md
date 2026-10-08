# Four-row affine products cannot create a half-symmetric profile from full cuts

This note checks a restricted but exact version of affine extraction.  Start
with a full fair cut system on `M>=4` source rows: every nonempty source cut (up to complement) is
present, with an arbitrary positive prime-power weight.  For the arithmetic
statement below, assume the Gaussian blocks for distinct cut classes have
disjoint rational-prime supports; shared-prime full-fair systems are outside
the scope of this note.  Form four output
rows by integer multiplicative affine combinations of the source rows.  The
rows of the coefficient matrix have sum one, so every output still has the
source radius before common-factor normalization.  The output is required to
have a fixed half-symmetric profile with two rows in each half.

The conclusion is a concrete no-go: four distinct output directions cannot
result.  This goes beyond deleting source rows or regrouping identical
restricted columns.  It uses the singleton cuts present in the full fair
matrix (and the pair cut containing any two singleton labels).  It does not
address balanced-only source matrices with more than four rows.

## 1. Exact affine valuation and phase formulas

Write a source block for a cut `T` as `kappa_T`, and orient its conjugate
pair by the indicator `1_T(i)`.  Let `A` be a `4` by `M` integer matrix with
each row sum one.  For the output row `r`, the exponent of `kappa_T` in the
formal product is

```text
n_r(T) = sum_i A[r,i] 1_T(i).                            (1)
```

The conjugate exponent is `1-n_r(T)`.  If the source block has Gaussian
prime exponent `e_T`, and distinct cut blocks have disjoint rational-prime
supports, then after clearing denominators and dividing the common Gaussian
gcd, its exact residual exponents are

```text
pi_T:     e_T (n_r(T)-min_s n_s(T)),
bar pi_T: e_T (max_s n_s(T)-n_r(T)).                    (2)
```

Thus the block contributes `e_T range_r(n_r(T)) log(p_T)/2` to the new
logarithmic radius.  This is the precise common-gcd cost in the
disjoint-support model; shared-prime layers can merge or cancel before a
min/max operation, so (2) is not asserted for them.

If `theta_i` are lifted source row arguments in an interval of width
`Delta`, the output phase difference is exactly

```text
theta'_r-theta'_s = sum_i (A[r,i]-A[s,i]) theta_i       (3)
```

modulo `2 pi`.  Complementary source cuts cause no ambiguity: because each
row of `A` has degree one, `n_r(T^c)=1-n_r(T)`, so the residual range and
the half-symmetric threshold condition are unchanged.  Since the coefficient difference has sum zero, its size on
the chosen lifts is at most
`||A[r,*]-A[s,*]||_1 Delta/2`.  Formulae (2)--(3) are the radius and angle
bookkeeping needed for any proposed normalized-length gain.

## 2. The four-row combinatorial obstruction

Fix the output partition `{1,2}|{3,4}`.  A nonconstant threshold layer in a
half-symmetric four-row profile has either one selected row in each half or
both rows in one half.  In other words, its support is a two-subset of the
four output positions.  Consequently, for every source cut `T`, the vector
`n(T)` must have its minimum and maximum each attained at least twice.  With
four coordinates this means either `n(T)` is constant or its values have
the form

```text
(u,u,v,v) after a permutation, with u != v.             (4)
```

The singleton source cuts force each column `A_{*,i}` to satisfy (4).  A
nonconstant column therefore determines one of the three pair partitions of
the four output positions.

Two nonconstant columns from different pair partitions cannot satisfy the
same condition after their sum.  After relabeling they have the form

```text
(u,u,v,v),       (r,s,r,s),       u != v, r != s.        (5)
```

Their sum is

```text
(u+r, u+s, v+r, v+s).                                  (6)
```

If `u-v` and `r-s` have the same sign, the first and fourth entries are the
unique maximum and minimum.  If their signs differ, the second and third
entries are the unique extrema.  Equality of either pair of candidate
extrema would force `u=v` or `r=s`; therefore the sum is never constant or
two-plus/two-minus.  The two-element source cut containing these singleton
indices then creates a forbidden singleton threshold layer.

It follows that all nonconstant columns of `A` have one common pair
partition. Every subset sum `n(T)` is then constant on the two pairs, so the
two output rows in each pair have identical coefficient rows and identical
valuations for every source block. The formal products themselves are equal,
even if the source rows carry arbitrary individual Gaussian units; no unit
convention is needed for this duplicate. Inserting extra row units after the
products would define a different model and is outside this conclusion.

Thus an integer affine product map from the full fair cut matrix cannot
produce four distinct rows satisfying the half-symmetric `k=2` profile.
The obstruction is exact for arbitrary positive block weights and prime
exponents. It does not assert that a balanced-only source matrix has the
same property.

## 3. Relation to larger output maps

The valuation range in (2) is the relevant cost for four-, six-, or
eight-row affine maps. For a balanced source benchmark with `B` blocks of
common norm weight `w`, the exact projective denominator is

```text
log H = (1/2) sum_T w_T range_r(n_r(T)),                (7)
```

where `w_T=e_T log(p_T)=log Norm(kappa_T)`.  This is exactly the
normalization in (2) of
[monomial_vector_subspace_height_barrier.md](monomial_vector_subspace_height_barrier.md):
there `E_p(v_r)=e_T n_r(T)`, so projective height is one half of the
primewise max-minus-min range.  Setting every `w_T=w` recovers the
equal-weight specialization.

The existing range theorem in
[monomial_vector_subspace_height_barrier.md](monomial_vector_subspace_height_barrier.md)
gives, for distinct projective characters, average range at least

```text
M=4:  q=3 -> 1/2,
M=6:  q=5 -> 3/5,
M=8:  q=7 -> 9/14,                                    (8)
```

where `q+1` is the number of output characters.  These are exact common
denominator costs after gcd normalization, not fair-norm assumptions about
an endpoint tuple.  They do not by themselves impose half symmetry; the
four-row argument above supplies the separate profile obstruction for a
full fair source.

Consequently, affine products/quotients offer no automatic extraction of
the half-symmetric theorem. Any successful six- or eight-row version would
need a new combinatorial argument controlling all source subset sums, plus
the exact radius and phase costs in (2)--(3). The present result makes no
uniform point-count claim and constructs no endpoint counterexample.

## Verification

The standard-library checker
[check_affine_full_fair_to_half_symmetric_obstruction.py](check_affine_full_fair_to_half_symmetric_obstruction.py)
enumerates all coefficient vectors in a bounded box, verifies the pair
partition sum obstruction exactly, and checks that the singleton and pair
source cuts force duplicate output classes.
