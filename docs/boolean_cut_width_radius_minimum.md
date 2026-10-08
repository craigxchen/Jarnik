# A sharp Boolean-width minimum for integral monomial row maps

For two-level prime allocations with equally weighted cuts, selecting k
source points minimizes the primitive squared radius among all integral
affine monomial maps with k distinct exponent rows. In particular, a
singular map that retains the full number of points cannot lower the
radius. This extends the nonsingular
[valuation-content floor](monomial_basis_valuation_content_floor.md),
but requires the additional equal-cut hypothesis.

Exact equal weights describe the limiting profile; Section 3 gives the
nonvacuous inequality for actual configurations with approximately equal
weights. Distinct prime supports cannot have exactly equal positive
logarithmic weights.

The result does not apply to general nested prime layers. Their images
can cancel inside one prime, as the existing
[nested reflection example](six_point_nested_reflection_budget.md)
demonstrates. It is not a uniform short-arc bound.

## 1. The discrete width theorem

Let a_1,...,a_k be distinct integer vectors in Z^n, all with the same
coordinate sum. For independent uniform bits X_j, put

```text
F(a)=E[max_i a_i.X-min_i a_i.X].
```

Then

```text
F(a)>=1-2^(1-k).                                      (1)
```

If the rows are affinely dependent, the stronger statement is

```text
F(a)>=1.                                              (2)
```

For k>=4, every configuration with F(a)<1 is, after a common integral
translation, a subset of the coordinate vectors or their negatives.
In this case equality holds in (1). Thus the gap (2) is strict relative
to the coordinate-simplex minimum. Equality in (2) occurs, for example,
at the four rows corresponding to the sets 12,14,23,34.

### Reduction to equal-size subsets

Write epsilon=2X-1. For a difference d=a_i-a_j, its coordinate sum
is zero, so

```text
F(a)>=E|d.X|=(1/2)E|d.epsilon|>=|d_j|/2.              (3)
```

The last step conditions on all signs except epsilon_j and uses
`(|t+d_j|+|t-d_j|)/2=max(|t|,|d_j|)`.
Consequently F(a)<1 forces every pairwise coordinate difference to
belong to {-1,0,1}. Subtracting the minimum of each column turns the
rows into distinct zero-one vectors of one common weight q. Regard
them as q-element subsets of [n].

For two such subsets at Johnson distance s (s elements removed and
s added), their expected projection distance is

```text
f_s=E|a_i.X-a_j.X|=s binom(2s,s)/4^s.
f_1=1/2, f_2=3/4, f_3=15/16, f_4=35/32.               (4)
```

This is the central binomial identity for the mean absolute value
of a 2s-step symmetric random walk, divided by two. It also follows
by summing `|r-s| binom(2s,r)/4^s`; the positive half telescopes
using `(r-s)binom(2s,r)=s[binom(2s-1,r-1)-binom(2s-1,r)]`.
The ratio `f_(s+1)/f_s=(2s+1)/(2s)` shows monotonicity.

For three rows there is the exact identity

```text
F(a,b,c)=(f_(s_ab)+f_(s_ac)+f_(s_bc))/2,               (5)
```

since the range of three real numbers is half the sum of their
three pairwise distances.

### Four rows force a coordinate simplex below width one

Suppose k>=4 and F(a)<1. Formula (4) excludes any pair at distance
at least four. A pair at distance three and any third row have
other distances at least one each and sum at least three. Equations
(4)--(5) then give width at least

```text
(15/16+1/2+3/4)/2=35/32>1.
```

Thus all pair distances are at most two. If some pair S,T has
distance two, every other row U must have distance one from both:
otherwise its triple has width at least
`(3/4+3/4+1/2)/2=1`.

Write `S=C union {1,2}` and `T=C union {3,4}`. Their common
distance-one neighbors are exactly `C union {a,b}`, where
`a in {1,2}` and `b in {3,4}`. Choose two distinct such neighbors.
After relabeling, the four rows, with their common C removed, are
one of

```text
{12,34,13,14},       {12,34,13,24}.                    (6)
```

Both have expected width exactly one. Here is an elementary count
on their four active bits: empty and full subsets contribute zero;
the eight subsets of size one or three contribute width one each;
the sum of widths over the six subsets of size two is eight. In the
first model those widths are 2,2,1,1,1,1; in the second they are
2,2,2,2,0,0. The total is sixteen out of sixteen choices. This
contradicts F(a)<1.

Every pair must therefore have distance one. Such a clique of
equal-size subsets is a star or a costar: either all contain a fixed
(q-1)-element subset, or all lie in a fixed (q+1)-element subset.
To verify this, fix adjacent rows C+a and C+b. Any other common
neighbor is either C+c or (C-d)+a+b. A new row of each type would
have distance two from the other. Hence only one type can occur.
The rows are consequently common translates of distinct e_j or
distinct -e_j. Their projection range is one except when all k
active bits agree, an event of probability 2^(1-k). This proves
the classification and (1) for k>=4.

### The small cases and affine dependence

For k=1, (1) is zero. For k=2, (3) and integrality give F>=1/2.
For k=3, each pair has expected distance at least 1/2 by the same
argument, so the three-number range identity gives F>=3/4.

No three distinct zero-one vectors are collinear: in any coordinate
where two differ, the line parameter must be zero or one for a
third zero-one vector. Therefore, when F<1, the binary reduction
already proves affine independence for k<=3. The classification
above proves it for k>=4. Its contrapositive is (2).

## 2. Exact consequence for two-level equal-cut configurations

Let m source Gaussian integers be primitive of common norm N.
Assume each split-prime column has the form `e_p=s_p 1_(S_p)` after
choosing an orientation. Let W_S be the sum of `s_p log p` over the
primes with unordered allocation cut S. Suppose each of the
`Q=2^(m-1)-1` nontrivial unordered cuts has the same weight w>0.
Then `log N=Qw`.

For an integral k-by-m row matrix A with A1=1, clearing denominators
and removing the exact Gaussian gcd gives

```text
log N_A=sum_S W_S width(A1_S)
       =w 2^(m-1) F(A).                               (7)
```

The second identity counts complementary cuts twice in the full
Boolean cube; the empty and full cuts contribute zero. Distinct
outputs require distinct exponent rows, since equal rows apply the
same monomial to the same source. Under the positive full-cut
hypothesis the converse also holds: a nonzero row difference d,
whose coordinate sum is zero, is detected by a singleton cut where
d_j is nonzero. Therefore (1) gives

```text
log N_A >= w(2^(m-1)-2^(m-k)).                         (8)
```

For k<=m, selecting k original rows attains (8), exactly as in the
[restriction gcd calculation](reflection_minimality_extraction.md).
For 4<=k<=m, equality in (8) requires the exponent rows to be common
translates of distinct coordinate vectors or their negatives. Such
common translation only multiplies all outputs by one rational
Gaussian factor; the negative case is common-norm conjugation
after a common factor. Primitive normalization removes that factor.

If the k rows of A are affinely dependent, (2) instead yields

```text
log N_A>=w 2^(m-1)=log N+w.                            (9)
```

In particular, any singular m-by-m affine row map with distinct
outputs increases the squared radius by at least exp(w). Any
nonsingular such map has N_A>=N, also following without equal-cut
weights from the valuation-content floor. Thus no integral affine
monomial map retaining all m distinct points lowers the radius
under the hypotheses of this section, regardless of entry size.

## 3. Stability and the nested-prime limitation

If `(1-epsilon)w<=W_S<=(1+epsilon)w`, the same singular/dependent
calculation gives

```text
log(N_A/N)>=w[1-epsilon(2^m-1)].                       (10)
```

Thus, for fixed m and sufficiently small relative cut error, the
strict obstruction survives even for arbitrarily large matrix
entries. Formula (10) still requires each prime column to have
only two levels.

For a general prime with nested layers `e=sum_t g_t 1_(S_t)` and
positive layer increments `g_t>0`, one
only has the upper triangle inequality

```text
width(Ae)<=sum_t g_t width(A1_(S_t)).
```

Equality need not hold: different layers at the same prime can
cancel after the row map. Replacing these layers by independent
binary primes would change the squared radius after transformation.
The exact nested countermodels therefore remain outside (7)--(10).

The [singular circuit tradeoff](singular_monomial_circuit_radius_tradeoff.md)
supplies a separate obstruction for bounded-height singular maps
using actual short-arc geometry, with no binary or fairness
assumption. Neither theorem supplies a point-count bound for an
arbitrary original cluster.

The proof is a prose theorem. An independent
[exact search and local audit](luna_boolean_width_search.md)
checks the two four-row patterns and the stated minima in bounded
integer boxes. These computations support the proof; they are not
used in place of its argument for arbitrary dimension.
