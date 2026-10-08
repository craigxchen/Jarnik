# A uniform nullity bound for squarefree cyclotomic indices

**Subsumed by the [all-odd-index theorem](cyclotomic_uniform_rank.md).**
The squarefree proof below was the intermediate step; prime-power
compression now removes its index restriction with the same constant.

There is a degree-independent bound for the one-class cyclotomic model
when every active index is squarefree, with no bound on its number of
prime divisors. The constant below is deliberately coarse.

Let `h` be a positive odd integer and `E` a finite set of odd squarefree
positive integers. Define

```text
M[d,e]=mu(e/d) if d|e, and 0 otherwise,    d odd, d<h,
D(E)=2*[1 in E]+sum_(e in E, e>1) phi(e).
```

Then

```text
D(E)<=4h  implies  dim_Q ker M<=1125.                  (1)
```

The proposed sharper bound of three is false, including in this
squarefree setting: see the
[exact nullity-four counterexamples](cyclotomic_nullity_four_counterexample.md).
The proof of (1) uses squarefree compression and a bounded-degree
intersection graph. The linked general theorem extends it to arbitrary
nonsquarefree active indices.
The empty-support case of (1) is immediate.

## 1. Compression to a divisor downset

Any row whose index is not squarefree is zero. On the other rows,
multiply row `d` and column `e` by their Möbius signs. The matrix becomes

```text
Z[d,e]=1 if d|e, and 0 otherwise.
```

This changes neither rank nor nullity. Fix a prime `p` and write
squarefree indices in lower/upper pairs `s,ps`, with `p` not dividing
`s`. Compress the column support by replacing every upper-only `ps`
by `s`; keep both members of an existing pair. The number of columns
is unchanged.

Split the evaluation rows into those not containing `p` and those
containing it. A lower and upper column over the same `s` have the
forms `(u_s,0)` and `(u_s,v_s)`. The old and compressed column spans
project onto the same span of the `u_s`. The kernel of this projection
on the old span contains every `(0,v_s)` arising from a retained
pair. The compressed span is precisely the direct sum of the span of
all `(u_s,0)` and the span of those paired vertical vectors. Hence

```text
rank(compressed matrix)<=rank(original matrix).
```

Thus nullity cannot decrease. The cost cannot increase either:
replacing `ps` by `s` removes the factor `p-1` from its totient when
`s>1`; for `s=1` it changes cost `p-1>=2` to `2`.

Apply such moves until none remains. Every nontrivial move decreases
the sum of the column indices, so the process terminates. Its final
support `E'` is divisor closed and satisfies

```text
D(E')<=D(E),
dim ker M_E <= dim ker M_E'.                           (2)
```

For a divisor-closed support, every row outside the support is zero.
The submatrix on its indices below `h` is upper triangular, ordered
increasingly, with diagonal entries one. Therefore

```text
dim ker M_E' = #{n in E':n>=h}.                        (3)
```

It remains to bound the count on the right.

## 2. A divisor-downset bound without a squarefree assumption

The following intermediate lemma permits arbitrary odd indices.
Let `E'` be a finite nonempty divisor downset of odd integers, and set

```text
U=sum_(e in E') phi(e),
r=#{n in E':n>=h}.
```

If `U<=4h`, then `r<=1125`.

To prove this, let `S` be the union of all complex roots of unity whose
orders lie in `E'`. Its cardinality is exactly `U`. For each `n in E'`,
the subgroup `H_n` consisting of the `n`th roots of unity lies in `S`,
since all divisors of `n` lie in the downset. In particular

```text
n<=U<=4h,
|H_n intersection H_m|=gcd(n,m).                      (4)
```

Make a graph on the `r` high indices, joining two distinct indices
when

```text
gcd(n,m)>=2h/15.
```

For a fixed `n`, write `g=gcd(n,m)`, `n=ag`, `m=bg`. The integers
`a,b` are positive, odd and coprime. By (4) and the edge condition,
both are at most `30`. There are only fifteen positive odd integers
at most `30`, so at most `225` possible pairs `(a,b)`, including
`(1,1)` for `m=n`. For fixed `n`, each pair determines `g=n/a` and
then `m=bn/a` uniquely. The graph consequently has maximum degree at
most `224`.

If `r>=1126`, a greedy independent-set algorithm finds six vertices:
each selected vertex removes at most `225` vertices, and five such
steps remove at most `1125`. Call the six independent orders
`n_1,...,n_6`. Every pair has `gcd(n_i,n_j)<2h/15`. The first two
terms of inclusion-exclusion give

```text
|union_(i=1)^6 H_(n_i)|
 >=sum_i n_i-sum_(i<j) gcd(n_i,n_j)
 >6h-15*(2h/15)
 =4h.
```

This contradicts containment in `S`, whose size is `U<=4h`.
The lemma follows. Applying it to (2)--(3), with
`U<=D(E')<=4h`, proves (1).

## 3. Consequences for the one-class polynomial model

Consider a family of distinct orientation polynomials with common
actual widths `W_e`, even base-layer counts, and pairwise contact
orders at least an odd `h`. Let

```text
L=sum_e phi(e) W_e<=4h,
E={e:W_e>0},
rho=dim ker M.
```

Assume every active `e` is squarefree. Multiplicities `W_e` remain
arbitrary: the squarefree assumption concerns the index, not its width.
The base width is at least two when active, so `D(E)<=L`. The theorem
therefore gives `rho<=1125`.

The [affine/obtuse bound](cyclotomic_affine_rank_reduction.md#5-an-affine-dimension-bound-for-every-whole-polynomial-fibre)
gives a whole-fiber count, without an arc-position assumption:

```text
number of polynomials <=2rho+1<=2251.                  (5)
```

If `L<4h`, the same proof improves this to

```text
number of polynomials <=rho+1<=1126.                 (6)
```

Indeed its affine embedding has pairwise inner products at most
`L-4h<0`. Such vectors are affinely independent: a nonzero affine
relation can be split into positive and negative parts with disjoint
supports; the squared norm of their common sum would then be a
strictly negative sum of cross inner products. Their affine dimension
is at most `rho`, proving (6).

For an actual growing fixed template of one resonant Fibonacci class,
the exact phase and allocation map in
[the affine bridge](cyclotomic_affine_rank_reduction.md#3-low-jets-and-the-actual-small-arc-count)
also gives at most `1126` distinct eventual points when its primitive
points lie on an arc of length at most `sqrt(2)*sqrt(R)`.
For general fixed `C>0`, subdivision gives the additional bound

```text
number of points <=1126 max(1,ceil(C/sqrt(2))).         (7)
```

The whole-fiber bound (5) is separately available under the same
one-class endpoint hypotheses, so (7) need not be used when it is
larger.

These bounds are independent of the layer degrees, rates, number of
distinct prime divisors of an index, and arbitrary layer
multiplicities, within the stated squarefree-index model. The fixed
template's asymptotic entry threshold may depend on the template.
No reduction of arbitrary circle configurations to this model is
asserted. Mixed-prime nonsquarefree indices are covered by the linked
general theorem; nonproportional affine classes remain outside it.
