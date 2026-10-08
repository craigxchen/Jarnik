# Complete balanced cuts block exact phase extraction

Consider extracting a smaller sign model by retaining selected rows and
blocks, and absorbing all omitted blocks into one common Gaussian
multiplier and the allowed Gaussian row units. This note gives the exact
criterion for that operation and a complete balanced-cut obstruction,
including regrouping identical or complementary restricted columns.
Other transformations or additional fixed algebraic row factors are
outside this criterion.

The obstruction is combinatorial. It is a statement about the full sign
pattern and arbitrary disjoint-prime block values; it is not an endpoint
configuration and does not prove that the balanced pattern is realized in a
short arc.

## 1. Exact regrouping of omitted blocks

Let `I` be a nonempty set of selected row labels and `J` a set of retained
block columns. For a sign matrix `S`, write

```text
Q_i = product_(j not in J) kappa_j^((1+S_ij)/2)
                         bar(kappa_j)^((1-S_ij)/2),
P_i = product_(j in J) kappa_j^((1+S_ij)/2)
                         bar(kappa_j)^((1-S_ij)/2).
```

For the disjoint-support conclusions below, assume every original block is a
nonunit (a varying unit column has no phase cost and is an exception).

Then the selected rows have the exact form `z_i=d u_i Q_i P_i`. The
omitted factors can be absorbed into one common Gaussian multiplier and
Gaussian row units precisely when `Q_i/Q_i0` is a Gaussian unit for every
`i in I`.
In that case `Q_i0` is absorbed into `d` and the remaining unit into `u_i`.

For every rational split prime choose one Gaussian prime `pi_p` and write
the conjugate-primitive block valuations as signed integers `e_jp`. The
omitted quotient has valuation

```text
1/2 sum_(j not in J) (S_ij-S_i0,j) e_jp              (1)
```

at `pi_p`, with the opposite value at `bar(pi_p)`. Consequently

```text
Q_i/Q_i0 is a unit for every i in I
iff sum_(j not in J) (S_ij-S_i0,j)e_jp = 0
   for every i in I and every p.                      (2)
```

The condition is exact for the actual blocks, including shared primes. Under
pairwise disjoint rational prime support, (2) is equivalent to

```text
S_ij = S_i0,j for every i in I and every j not in J.  (3)
```

Thus an omitted nonconstant column really is a moving row phase. It cannot
be discarded by appealing to a row unit: a conjugate-primitive nonunit is
not associate to its conjugate.

## 2. Complete balanced cuts

Put `M=2k>=6`, and index columns by all `k`-subsets `T` of `[M]`:

```text
S_i(T) = +1 if i in T, and -1 otherwise.
```

Let `D` be the row-difference matrix. Every column of `S` sums to zero, so
`Dq=0` implies `Sq=c 1` and hence `c=0`. Also

```text
S_0 = -(1/M) sum_(i=1)^(M-1) (S_i-S_0).
```

Therefore

```text
ker(D)=ker(S),       V:=rowspan(D)=rowspan(S).          (4)
```

The full kernel has no signed coordinate group of the type needed by the
phase-group theorem. More strongly:

> **Two-support lemma.** No nonzero vector supported on at most two columns
> belongs to `V`.

To prove this, suppose `x in V` is supported on at most two columns. A
row-space representation gives `f(T)=sum_i a_i S_i(T)`, with `f(T)=0`
outside those two columns. Fix `i!=j`. There are
`binom(M-2,k-1)>=6` sets `C`, of size `k-1`, disjoint from `{i,j}`. At most
two of the corresponding pairs `C union {i}`, `C union {j}` can be
exceptional, so choose `C` for which both values of `f` vanish. Their
difference is `2(a_i-a_j)`, hence `a_i=a_j`. All `a_i` are equal, and then
`f(T)=0` for every balanced set because `sum_i S_i(T)=0`. This gives `x=0`,
a contradiction.

In particular, no singleton `e_j` lies in `V`, and no signed pair
`epsilon_j e_j-epsilon_h e_h` lies in `V`. A group condition

```text
r_I = sum_(j in I) epsilon_j e_j in ker(D),
epsilon_j e_j-epsilon_h e_h in V                  (5)
```

therefore has no nonempty solution: groups of size at least two violate the
two-support lemma, while a singleton would require `e_j in ker(D)`, and the
`j`-th sign column is nonzero. This rules out the full-space signed partition
extraction before any Gaussian arithmetic or arc estimate is used.

## 3. Regrouping identical restricted columns does not rescue extraction

Beyond deleting original columns, one may regroup columns whose
restrictions to `I` agree up to complement. For
one such class, orient the blocks consistently and multiply them into one
effective block. With disjoint rational prime supports this effective block
is still conjugate-primitive and is a nonunit whenever the class is
nonempty. Constant restricted classes are absorbed into the common factor.

The number of nonconstant classes after this optimal regrouping is

```text
b_red(M,t) = (1/2) sum_(h=max(1,t-k))^(min(t-1,k)) binom(t,h).  (6)
```

Here `h` is the number of selected rows carrying the plus sign. A completion
to a balanced `k`-set exists exactly when
`max(0,t-k)<=h<=min(t,k)`, and quotienting by complement identifies `h` with
`t-h`. Thus (6) counts precisely the effective nonunit blocks, rather than
the original columns.

The first values are `b_red(M,1)=0`, `b_red(M,2)=1`, and
`b_red(M,3)=3`. For `t=4` all seven nonconstant classes occur. For even
`t>=6`, the middle slice alone gives
`b_red(M,t)>=binom(t,t/2)/2>t`; for odd `t>=5`, the two central slices give
`b_red(M,t)>=binom(t,(t-1)/2)>t`. The central slice is feasible even when
`t>k`, since `t<=2k`. Consequently a rank-one `t'`-column profile would
require `t'-1<=t-1`, hence `t'<=t`, while regrouping leaves `t'>t` for every
`t>=4`; `t<=3` cannot provide nine blocks.

There is a parallel post-regrouping group obstruction. On odd `t>=5`, use
the slice `h=(t-1)/2`. A quotient class has one representative there, and a
two-support row-space vector can mark at most two representatives. For any
two row labels, there are
`binom(t-2,h-1)>2` exchanges `C union {i}`, `C union {j}`; one avoids both
exceptions and forces the two coefficients to agree. On even `t>=6`, use
the middle slice `h=t/2`. Each quotient exception has two complementary
representatives, so at most four slice sets are exceptional, while
`binom(t-2,t/2-1)>=6`; the same exchange argument applies. Here the
row-difference space is used: write a certificate as
`f(T)=sum_i a_i S_i(T)` with `sum_i a_i=0`. The case `t=4` is the finite
seven-class system. Its singleton classes have values `2a_i`; if a nonzero
`f` had at most two nonzero classes, the zero singleton classes would leave
at most two nonzero coefficients. Since `sum_i a_i=0`, there must be
exactly two, with opposite values. Their mutual pair class is zero, but
the other two pair classes are nonzero. Together with the two nonzero
singleton classes this gives four nonzero classes, a contradiction.
Thus no nonzero row-difference
vector is supported on at most two reduced classes for `t>=4`. The small case
`t=3` really does have a three-column rank-one profile, but it has too few
blocks for the nine-block Roth theorem.

This proves an exact no-go for extracting a rank-one phase profile with at
least nine effective nonunit blocks from a complete balanced-cut matrix by
selecting rows, regrouping restricted columns, and absorbing the remaining
columns. It also handles the natural signed-group generalization after
regrouping.

## 4. Scope of the obstruction

The disjoint-support hypothesis, together with nonunit original blocks,
makes every nonconstant regrouped class a nonunit conjugate-primitive
effective block and gives the equivalence between (2) and the columnwise
condition (3). With shared primes, omitted columns can cancel in (2), and
that cancellation must be checked prime by prime; it cannot be inferred
from the sign matrix. Conversely, (2) itself remains an exact criterion for
the chosen actual blocks.

The balanced profile also has fractional linear-phase packing optimum
`2(M-1)/(M-2)<4`, as recorded in
[`coupled_block_roth_packing.md`](coupled_block_roth_packing.md). The present
argument is a different obstruction: it blocks the hoped-for bounded
submatrix extraction before packing is considered. Neither fact is a
counterexample to the endpoint theorem. Fixed-pattern exclusion, even with
an effective bound for each fixed pattern, still does not give a uniform
point-count bound as the pattern varies.

## Verification

The standard-library checker
[`check_complete_balanced_omitted_column_obstruction.py`](check_complete_balanced_omitted_column_obstruction.py)
enumerates the balanced columns for `M=6,8`, verifies (4), checks every
one- and two-support vector against the exact row space, and checks both the
original and regrouped counting formulas. It also checks the reduced-pattern
support obstruction. It only certifies finite linear algebra; the regrouping
identity and the disjoint-prime implication are exact symbolic arguments.
