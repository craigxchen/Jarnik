# A sharp sign-variation test for the degree-22 odd-column code

This note records one bounded consequence of the nonlinear odd-moment
equations for the code in
`tensor24_odd_column_certificate_gap_fixture.json`.  It also gives an exact
obstruction to using that consequence alone to exclude the code.  It does not
construct an odd-moment solution.

Let $S$ be the stored $8\mathbin{\times}22$ sign matrix and put

\[
 B_i=(S_i-S_7)/2,\qquad 0\leq i<7.
\]

Suppose that nonzero real numbers $a_j$, with pairwise distinct absolute
values, obey

\[
 S(a^{2k+1})=\alpha_k\mathbf 1,\qquad 0\leq k\leq4.
\]

Then $B(a^{2k+1})=0$.  For every nonzero vector $c$ in the real row
space of $B$, set $t_j=a_j^2$ and $w_j=c_ja_j$.  Ordered by increasing
$t_j$, the nonzero sequence $w_j$ must have at least five sign changes.
Indeed,

\[
 \sum_j w_jt_j^k=0,\qquad 0\leq k\leq4.
\]

If the sequence had at most four sign changes, choose one point strictly
between consecutive $t_j$ values at each change.  Up to an overall sign, the
product of the corresponding linear factors is a polynomial $q$ of degree
at most four with $w_jq(t_j)>0$ at every nonzero $w_j$.  This contradicts
$\sum_jw_jq(t_j)=0$.

For this code the condition is sharp but compatible.  With zero-based column
indices, take the increasing-$t$ order

```text
6, 17, 8, 14, 4, 9, 18, 21, 0, 7, 12,
3, 11, 15, 16, 10, 2, 1, 5, 13, 20, 19
```

and the column-indexed signs $\operatorname{sgn}(a_j)$

```text
+ + + + + - + + + - - - + + - + - + + - - +
```

Every nonzero vector in the row space then has at least five sign changes.
It is enough to check the support-minimal row-space vectors.  Every nonzero
row-space vector contains a conformal support-minimal vector.  To see this,
intersect the row space with the cone of vectors having its weak coordinate
signs, and set its zero coordinates to zero.  In that cone, move toward a
boundary until another coordinate vanishes, continuing to an extreme ray.
The cone is pointed because it lies in one coordinate orthant; an extreme
ray has minimal support, since a smaller supported vector would permit a
small positive and negative perturbation within the cone.
Between two consecutive nonzero entries of the smaller vector with opposite
signs, the larger vector must change sign at least once.  Hence adding its
other nonzero entries cannot reduce the number of changes.

The exact checker enumerates all support-minimal vectors without a coefficient
bound.  A support-minimal vector has a maximal proper zero-column flat.  Since
$B$ has rank seven, that flat has rank six.  Choosing six independent
columns in it determines a normal vector through its signed $6\mathbin{\times}6$
minors.  Conversely every nonzero normal obtained this way has precisely the
maximal zero flat belonging to a support-minimal vector.  Deduplication up to
overall sign gives 10,155 vectors, with support distribution

```text
support:   7   10   11   12   13   14    15    16
count:     1    1   44  105  154  639  2086  7125
```

Their minimum variation under the displayed order and signs is exactly five.
The unique support-seven vector, in stored column order, is

```text
(0,0,0,0,0,0,0,0,0,0,3,-1,0,0,0,-2,-3,-1,0,-2,0,-2).
```

Thus positivity of the squared nodes gives a real simultaneous Chebyshev
sign-variation constraint, but the oriented matroid of this code satisfies
it.  This does not rule out stronger moment-matrix arguments.  Any exclusion
of the full nonlinear system must use coefficient magnitudes or additional
polynomial relations, not sign variation alone.

The [critical root-order theorem](critical_odd_moment_root_order.md)
subsequently supplies such an additional polynomial constraint. Its
simultaneous prefix test [excludes this fixed code](tensor24_odd_column_critical_order_obstruction.md)
and [the entire normalized odd-positive tensor family](tensor24_odd_column_family_critical_order.md).
The compatible order recorded here satisfies only the weaker variation
condition.

Run

```text
python3 docs/check_tensor24_odd_column_variation.py
```

for the integer-minor enumeration and the displayed compatibility check.
