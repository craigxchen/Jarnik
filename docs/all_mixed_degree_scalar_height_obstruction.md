# Every degree-one/two factor has the scalar height obstruction

Changing the multidegree does not produce a strict successive-minima
gain from the scalar balanced-cut congruences. This holds for **every**
positive degree vector with entries in `{1,2}`, in every row count, and
does not assume that its weighted-balanced evaluation map has full rank.
These are exactly the supported multidegrees that can arise by factoring
a degree-two-each relation in
[minimal_support_invariant_descent.md](minimal_support_invariant_descent.md).

This is an audited obstruction to a specified method, not an exclusion of
the extracted near-uniform profile. It does not give a lower height for
an individual irreducible zero in the balanced kernel. In particular,
the five-row degree-two case has no weighted-balanced cuts at all.

The later [arbitrary-degree theorem](arbitrary_degree_scalar_height_obstruction.md)
extends this obstruction to all admissible positive integer multidegrees.
The specialized table and argument below remain valid.

## 1. The actual arithmetic inputs and the most favorable denominator

Use the actual full-cut Gaussian system from
[invariant_degree_variation.md](invariant_degree_variation.md), with a
fixed number `m` of supported rows, inherited block logarithms
`w+o(w)`, and correction and primitive pair-residue logarithms `o(w)`.
For positive degrees `d_i in {1,2}`, put `D=sum_i d_i` and assume `D`
is even. An invariant of odd total degree is zero. Write

```text
r = [x^(D/2)] product_i(1+x+...+x^d_i)
    -[x^(D/2-1)] product_i(1+x+...+x^d_i),
s = (1/2) #{S subset [m]: sum_(i in S)d_i=D/2},
A = D 2^(m-3)-sum_(S subset [m]) max(0,d(S)-D/2).       (1)
```

The primitive numerical evaluation vector has logarithmic height
`Aw+o(w)`. This is the exact graph-content formula already proved for
unequal degrees; it retains the inherited full-support block until its
common content is divided out. Complementary scalar congruences force
any numerical zero outside their kernel to have coefficient logarithm
at least `2w-o(w)`.

Let `rho` be the actual rank of the scalar map, and let `K` be its
kernel. Only `rho<=min(r,s)` will be used. Numerical evaluation `ev` is
a nonzero functional: a nonzero bracket monomial has nonzero value at
distinct directions. Its integer relation lattice has rank `r-1`.

If `ev|K` is nonzero, the relations in `K` have dimension `r-rho-1`.
The first successive minimum guaranteed to leave `K` is consequently
index `r-rho`, whose Minkowski-product denominator is `rho`.

If `ev|K=0`, that internal dimension is `r-rho`. The required index is
`r-rho+1`, with denominator `rho-1`; when `rho<=1` there need not be
any relation outside `K`. This includes `K=0`, where the denominator is
`r-1`.

Thus in every case where this argument supplies an outside relation its
denominator `b` satisfies

```text
1 <= b <= min(s,r-1).                                  (2)
```

The guarantee from the full relation covolume is `Aw/b+o(w)`. To beat
the complementary scalar threshold by this argument requires `A/b<2`.
Using the most favorable possible denominator in (2), even without
knowing the actual scalar rank, the next theorem rules this out.

## 2. A uniform combinatorial inequality

For every `m>=4` and every positive `d_i in {1,2}` of even total degree,

```text
A >= 2 min(s,r-1).                                     (3)
```

Whenever `s>0`, the inequality is strict for `m>=5`. At `m=4`, all three
possible unordered degree patterns have equality in (3):

```text
degrees       r   s   A   min(s,r-1)
1,1,1,1       2   3   2       1
1,1,2,2       2   2   2       1
2,2,2,2       3   3   4       2
```

There are no possible numerical zero factors on at most three distinct
rows, since each nonzero fixed-multidegree invariant there is a scalar
multiple of one bracket monomial.

Here is a proof for all large `m`, followed by a short exact certificate
covering the remaining finitely many degree patterns. Choose independent
uniform signs `epsilon_i`. Complement symmetry in (1) gives

```text
A = 2^(m-3) [D-2 E|sum_i d_i epsilon_i|].               (4)
```

Cauchy--Schwarz and `d_i^2<=2d_i` imply

```text
A >= 2^(m-3) [D-2 sqrt(2D)].                           (5)
```

Also `2s<=2^(m-1)`: pairing subsets that differ only in membership of
one fixed row, at most one member of a pair has weight `D/2`, because
that row has positive degree. If `m>=16`, then `D>=16`. The increasing
function `D-2sqrt(2D)` is strictly greater than four for `D>=16`:
at sixteen this is `16-8sqrt(2)>4`, since `sqrt(2)<3/2`.
Hence

```text
A > 2^(m-1) >= 2s,             m>=16.                  (6)
```

For `4<=m<=15`, a degree vector is determined up to permutation by its
even number `a` of ones, with `m-a` twos. The following table is the
exact minimum of `A/min(s,r-1)` among patterns with positive denominator.
It is computed from the finite binomial sums (1); the accompanying
standard-library checker exhausts **all** patterns in this range.

| m | Minimum | A pattern attaining it: number of ones, twos |
|---:|---:|---:|
| 4 | 2 | 0, 4 |
| 5 | 8/3 | 2, 3 |
| 6 | 18/5 | 0, 6 |
| 7 | 29/6 | 4, 3 |
| 8 | 232/35 | 0, 8 |
| 9 | 8 | 6, 3 |
| 10 | 209/21 | 6, 4 |
| 11 | 49/4 | 8, 3 |
| 12 | 4200/323 | 8, 4 |
| 13 | 4901/300 | 8, 5 |
| 14 | 20250/1229 | 10, 4 |
| 15 | 23100/1163 | 10, 5 |

Every entry after four is strictly greater than two, proving the
remaining cases of (3). These are the most favorable abstract
denominators, not assertions that the corresponding scalar maps attain
them. A smaller actual rank makes the guarantee weaker.

## 3. Exact consequence for irreducible-support descent

Suppose an ambient `M`-row relation has an irreducible numerical-zero
factor supported on `m` rows. Its coefficient logarithm increases by
only `O_M(1)`, while the inherited profile scale on those rows is

```text
w_m=2^(M-m) w.                                        (7)
```

This amplification is real and useful when there is an independent
lower coefficient bound for that individual factor. Nothing above
removes it. The theorem only says that **restarting** the full-covolume
successive-minima construction in any of the possible residual
multidegrees cannot itself establish a strict comparison with the
scalar threshold: both its guarantee and that threshold use `w_m`,
and (3) still applies.

Some unequal choices do improve the numerical guarantee over equal
degree two; for example the optimistic ratio at ten rows is `209/21`
instead of `1300/126`. They remain far above two. The obstruction is
therefore stronger than a claim about the equal-degree family, while
leaving all-cut image compatibility, individual irreducible height
bounds, and source-prime residue information open.

The near-uniform extraction by itself does not assert the negligible
correction and residue hypotheses used in Section 1 for every fixed
tuple. The result is conditional on that actual arithmetic model, and
shows that even this stronger favorable model supplies no strict gain
from the specified scalar/covolume package. It does not turn the
endpoint profile exclusion into a proved statement.
