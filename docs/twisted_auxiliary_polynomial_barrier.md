# Moving small twists and the auxiliary-polynomial comparison

This note extends the degree-versus-divisibility comparison in
`hermitian_polynomial_barrier.md` to the rounded endpoint-allocation
system. The extension permits moving coefficients, arbitrary bounded
total degree, and conjugation. It remains a restriction on a specified
proof method, not a theorem excluding all polynomial arguments.

## 1. The rounded system restores the archimedean inputs

Fix `m` and let the scale `W` tend to infinity. For every subset
`S` of `[m]`, let `H_S` be a Gaussian block. Assume distinct blocks
and their conjugates are pairwise coprime, and

```text
log N(H_S)=W/2^m+o(W).
```

Put

```text
A_i=product_(S containing i) H_S.
```

The rounding argument produces Gaussian integers `V_i!=0` such that

```text
P_i=V_i A_i=X_i+iY_i,
log max(1,|V_i|)=o(W),
Y_i!=0,       log max(1,|Y_i|)=o(W).
```

The original angular precision is retained: the arguments of `P_i`
lie within `O(exp(-W/4))` of a common real ray after choosing signs.
Thus

```text
log |P_i|=W/4+o(W),
log |X_i|=W/4+o(W).
```

The empty block is allowed. The scale here is the inherited original
`W`; changing radius after deleting the empty block must also change
the angular exponent, as in the preceding Hermitian note.

For `i in S`, both prescribed divisibilities remain true:

```text
H_S divides P_i,
conjugate(H_S) divides conjugate(P_i).
```

Small row twists therefore restore the small imaginary coordinates
needed by an auxiliary-polynomial argument without losing the leading
prescribed block divisibility. They do not make the rounded `A_i`
themselves nearly real.

## 2. A moving-coefficient proposition

Let `Q_W(X_1,...,X_m)` be a nonzero polynomial over `Q(i)`, of total
degree `d_W<=D`, where `D` is fixed. Suppose every coefficient has
absolute logarithmic height `o(W)`. There are only boundedly many
coefficients, so one common denominator of height `o(W)` makes the
polynomial Gaussian-integral. After this clearing, all coefficient
moduli are `exp(o(W))`.

For a subset `S`, define `nu_S^+` and `nu_S^-` to be the generic
orders along the coordinate subspaces

```text
X_j=-iY_j for j in S,
X_j=+iY_j for j in S,
```

respectively. The empty-set orders are zero. These are ordinary
characteristic-zero ideal-membership orders for each polynomial in
the sequence. They may vary with `W`.

Then

```text
sum_S (nu_S^+ + nu_S^-) <= d_W 2^(m-1).           (1)
```

Consequently the logarithmic modulus of the divisor obtained from
these generic prescribed block orders is at most

```text
(1/2) sum_S log N(H_S)(nu_S^+ + nu_S^-)
 <= d_W W/4+o(W).                                 (2)
```

The corresponding archimedean estimate is

```text
|Q_W(X_1,...,X_m)| <= exp(d_W W/4+o(W)).           (3)
```

Thus there is no fixed positive exponent margin between (2) and (3).
This is uniform over the stated class of moving polynomials; it is
not restricted to one fixed coefficient vector.

### Proof of the moving-grid inequality

Use the rectangular grid whose coordinate sets are
`{-iY_j,+iY_j}`. Each set has two distinct elements because `Y_j!=0`.
At the vertex with the negative choice on `S` and the positive choice
on `S^c`, the Taylor expansion has multiplicity at least
`nu_S^+ + nu_(S^c)^-`. Indeed the two prescribed orders involve
disjoint variable sets, so every Taylor monomial has at least their
sum of degrees.

The multiplicity grid inequality gives an upper bound
`d_W 2^(m-1)` for the sum of multiplicities. Its elementary proof by
induction on the number of variables works over `Q(i)` for any
two-element coordinate sets, including these moving sets: write the
polynomial in the last variable, use the multiplicity of its leading
coefficient at a smaller-grid point, and sum the remaining univariate
multiplicities, whose total is at most the last-variable degree.
Reindexing complements proves (1).

All orders are at most `D`, and there are only `2^m` subsets, so the
`o(W)` profile errors in (2) remain `o(W)` even when the polynomial
and its orders vary. Estimate (3) follows from bounded degree,
boundedly many monomials, the coefficient-height assumption, and
`|X_i|=exp(W/4+o(W))`.

## 3. Coefficient primes do not hide a leading generic gain

For a fixed subset and sign, expand the denominator-cleared `Q_W`
in the translated variables `X_j+iY_j` or `X_j-iY_j`, translating
the other variables as well if convenient. Every translated coefficient
is a bounded-degree polynomial in the original coefficients and the
`Y_j`. Hence it is a Gaussian integer of modulus `exp(o(W))`.

Choose one nonzero translated coefficient `b_(S,+)` whose total
degree in the variables indexed by `S` is exactly `nu_S^+`.
At a Gaussian prime `pi` occurring to exponent `e` in `H_S`, the
generic valuation supplied by the prescribed divisibilities is
the minimum of the coefficient valuation plus the weighted monomial
valuation. In particular it is no larger than

```text
e nu_S^+ + v_pi(b_(S,+)).                         (4)
```

The same statement holds at the conjugate place and the opposite
translated expansion. Summing the second terms in (4), over all
conductor primes and then the finitely many subsets and signs, costs
at most

```text
sum_(S,sign) log max(1,|b_(S,sign)|)=o(W).
```

Thus allowing the coefficients to be specially divisible at moving
conductor primes does not change the leading comparison. The
denominator-clearing factor also has logarithmic modulus `o(W)`.

The twists `V_i` can provide additional known divisibility beyond
that of `A_i`. A degree-`D` monomial uses at most `D` copies of any
such extra valuation. Their total available logarithmic modulus is
bounded by a fixed-degree multiple of
`sum_i log max(1,|V_i|)=o(W)`, with conjugate factors handled in the
same way. They cannot supply a fixed leading gain solely through
their own prime factors.

Here “generic” refers to the valuation obtained by expanding the
polynomial and using these prescribed divisibilities and coefficient
valuations. Special cancellations after substituting the actual
outside coordinates can make the value more divisible, or zero.
Those effects are not bounded by (4) and remain possible extra
arithmetic inputs.

## 4. Rank, Gram, and conjugation expressions fit the proposition

Every bounded-degree polynomial expression in

```text
P_i, conjugate(P_i), Y_i, V_i, conjugate(V_i)
```

with negligible-height coefficients becomes, after the exact
substitutions

```text
P_i=X_i+iY_i,
conjugate(P_i)=X_i-iY_i,
```

a polynomial in the real coordinates `X_i` with negligible-height
coefficients and bounded degree. Rational small-height data can be
cleared with the same bounded-degree height cost. If the resulting
polynomial is nonzero, Sections 2--3 apply.

For example, the Gram entries are

```text
P_i conjugate(P_j)
 = X_i X_j+Y_i Y_j+i(Y_i X_j-X_i Y_j).
```

The universal rank-one minors and cocycles become identically zero
after substitution. They therefore do not produce nonzero integer
values to which a product-formula contradiction could be applied.
Other expressions can be nonzero and fall under the comparison above.

This observation handles moving bounded-degree Hermitian expressions
more broadly than the fixed multihomogeneous presentation. In that
presentation, the Taylor order `v` simply reduces the maximal degree
in the large coordinates from `D` to `D-v`, because every `Y_i`
has become a small-height coefficient here.

## 5. Precise remaining possibilities

The rounded normal form does not invalidate the earlier polynomial
barrier: its small twists restore exactly the archimedean inputs, and
the moving-grid argument shows that their height cost does not create
a strict generic divisibility margin.

This result does not cover every use of the norm data. For example
`N(P_i)=N(V_i)N(A_i)` contains a conductor-sized quantity. Treating
that quantity as a polynomial coefficient violates the negligible-
height hypothesis unless its size and divisibility are tracked
separately. Likewise a cancellation specific to the actual block
values is stronger information than generic coordinate-subspace
membership. Either feature could be part of a different argument.

What has been established is an exact obstruction to deriving a
positive exponent saving solely from bounded-degree moving small
coefficients, the small twists and residues, generic prescribed
conductor divisibility, and the usual polynomial size estimate.
No new impossibility theorem for the rounded system itself is proved.

An exact Gaussian-coefficient calculation checked (1) on 400
inhomogeneous polynomials with moving nonzero centers, using two
through five variables and products of shifted isotropic factors,
norm factors, weighted brackets, and affine coordinate factors.
There were 217 equality cases and no violation. This finite check
supplements the moving-grid proof and is not a Lean formalization.
