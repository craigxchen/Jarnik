# Singular monomial maps pay for retaining distinct short-arc points

Singular affine monomial maps can erase prime-allocation columns, unlike
the nonsingular maps covered by
[the valuation-content floor](monomial_basis_valuation_content_floor.md).
They also force a multiplicative relation among their outputs. The
existing product-rigidity theorem converts that relation into an exact
radius-versus-angle obstruction.

In particular, using cross products across a prime cut to restore all k
source points cannot lower the radius once the source normalized arc
constant is at most `2/sqrt(k)`. The proof uses the actual angles and
cardinality, independently of all prime weights.

## 1. A relation gives a finite radius tradeoff

Let the primitive Gaussian source tuple have common norm N and angular
span Delta; write its normalized arc constant as

```text
C = Delta N^(1/4).
```

Choose an integer row matrix A with row sums one and raw outputs
`w_i=product_j z_j^(A_ij)`. Clear one common denominator and remove
the exact common Gaussian gcd, obtaining distinct Gaussian outputs
`w'_i` of common norm N'. Suppose their actual angular span is at
most `K Delta`, where K is a specified positive constant.

If an integer vector c satisfies

```text
c^T A=0,      c!=0,      h=sum_(c_i>0)c_i,
```

then `sum c_i=0` and the outputs obey an h-versus-h multiplicative
relation. Common normalization preserves it because the two sides
have equal total degree. Distinctness of the outputs makes the
multisets in that relation different. By
[angular moment product rigidity](angular_moment_product_rigidity.md),
their normalized arc constant C' must satisfy

```text
h (C')^2 > 8.
```

But `C'<=K C (N'/N)^(1/4)`. Therefore

```text
N'/N > 64 / (h^2 K^4 C^4).                              (1)
```

The strict inequality includes the equality case of the B_h theorem.
This argument does not need multiplicative independence of the source.
It does need distinct outputs; otherwise the displayed relation might
be merely an equality of multisets.

## 2. Every bounded-height singular matrix supplies such a relation

Suppose A is k-by-k, singular, and its rows have l1 norm at most L.
Take a minimal linearly dependent set of q rows, so `q<=k` and its
rank is `q-1`. Choose `q-1` columns of that rank. Their signed maximal
minors give a nonzero integer dependence c on those rows. Hadamard's
inequality bounds each minor by `L^(q-1)`. After dividing the common
integer content, this gives

```text
h=||c||_1/2 <= q L^(q-1)/2 <= k L^(k-1)/2.               (2)
```

The chosen columns span the column space of the q-row matrix, so
this dependence annihilates every column of A, not only the selected
ones. Its coefficient sum is zero because every row sum of A is one.

For these monomials, a source angular lift also gives

```text
Delta' <= K_A Delta,
K_A = (1/2) max_(i,j) ||A_i-A_j||_1 <= L.                (3)
```

To prove (3), subtract any two output arguments. Their coefficient
vector has sum zero; on an interval of length Delta, its absolute
value is bounded by half its l1 norm times Delta. The resulting
lifted arguments all lie in an interval of the stated length.

Combining (1)--(3), any such matrix retaining k distinct outputs has

```text
N'/N > 256 / (k^2 L^(2k-2) K_A^4 C^4)
     >= 256 / (k^2 L^(2k+2) C^4).                       (4)
```

Consequently a nonincreasing radius requires

```text
C > 4 / (sqrt(k) L^((k+1)/2)).                          (5)
```

For fixed k and bounded L, singular row changes therefore cannot
retain k distinct outputs, preserve a comparably short arc, and
decrease the radius along a family with C tending to zero. Matrices
whose entries grow with the configuration remain outside this
bounded-height conclusion.

If k-1 rows are retained and the missing row is reconstructed as
one affine monomial in them, there is a simpler bound. If that
monomial has positive coefficient mass P, its defining relation
already has length P. Thus distinct outputs require `P(C')^2>8`,
without the determinant estimate (2).

## 3. Cross products across a cut have a short cycle certificate

Partition the k source labels into nonempty parts S,T. Fix one
source label s and choose outputs from

```text
w_ij = z_i z_j/z_s,              i in S, j in T.         (6)
```

If a source prime has a binary allocation constant on S and on T,
all outputs in (6) have the same valuation at that prime. Primitive
normalization therefore erases that column. This is the intended
arithmetic advantage of crossing the cut.

Associate a selected output to the edge ij of a bipartite graph on
the k source labels. Any k selected edges contain a simple cycle:
a forest on k vertices has at most k-1 edges. Write the cycle length
as `2h<=k`. Its alternating edge products satisfy an exact
h-versus-h relation, because each source vertex and the common
denominator z_s occur equally often on both sides.

For k distinct selected outputs, the two alternating multisets are
different. The actual angular span of (6) is at most
`span(S)+span(T)<=2 Delta`, so (1) with `K=2` gives

```text
N'/N > 4/(h^2 C^4) >= 16/(k^2 C^4).                    (7)
```

In particular,

```text
C<=2/sqrt(k)  ==>  N'>N.                                (8)
```

If a chosen list contains duplicate outputs, it already fails to
restore k distinct points. Thus either cardinality fails or (7)
holds; there is no assumption that all cross products are distinct.
When the source is a B_2 set, distinct edges in this bipartite graph
do indeed give distinct products, but that extra observation is
unnecessary for the dichotomy.

One may also apply the argument to any selected subgraph containing
a short cycle, even if fewer than k outputs were selected. A
four-cycle gives h=2; stronger existing rectangle separation bounds
can improve the numerical constant in that special case.

## 4. Scope

These inequalities show why two natural singular repairs do not
complete a backwards descent: bounded-height rank reduction carries
a bounded multiplicative circuit, while restoring k outputs by
crossing a prime cut necessarily introduces a graph cycle. The
short-arc product theorem forces the resulting radius cost.

They do not exclude singular maps of unbounded coefficient height,
nonmonomial replacements, or a strategy that deliberately loses
points and gains enough radius for a separate induction. No uniform
point count or improvement in its radius dependence is proved here.
