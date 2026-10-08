# Why product rigidity does not settle the balanced-pattern obstruction

The recent product-rigidity theorems are valid restrictions on endpoint
arcs. This audit identifies a limitation of that route: the complete
balanced cut patterns already have no multiplicative relations of any
length. Consequently increasing the permitted relation length cannot
by itself exclude those patterns.

## Full affine rank

Let `M>=2` be even, `B=binom(M,M/2)`, and let `A` be the `M`-by-`B`
matrix whose columns are all indicator vectors of subsets of size
`M/2`. Each diagonal entry of `A A^T` is `B/2`. Each off-diagonal
entry is

```text
binom(M-2,M/2-2) = B(M-2)/(4(M-1)),
```

with value zero when `M=2`. Thus for every real vector `c` with
`sum c_i=0`,

```text
||A^T c||_2^2 = B M/(4(M-1)) ||c||_2^2.           (1)
```

In particular, `sum c_i=0` and `A^T c=0` imply `c=0`. The rows
are affinely independent. Every subset of `h` rows is affinely
independent too, by extending any purported dependence by zero.
This holds over the reals, so there are no integer circuits with
large coefficients either.

Subtracting a coordinate's minimum over a row subset does not change
an affine relation, since its coefficients sum to zero. Therefore
the same independence survives the common-factor normalization of
that subset. Restriction and primitive normalization cannot turn
eight of these rows into a width-four exponent family: eight
affinely independent rows require at least seven active coordinates.

## Realization by actual Gaussian integers

Assign to each column `S` a different rational prime `p_S=1 mod 4`
and choose a Gaussian prime `pi_S` of norm `p_S`. Define

```text
z_i = product_(S contains i) pi_S
        product_(S does not contain i) conjugate(pi_S).
```

All these Gaussian integers have squared modulus `product_S p_S`.
They have no common Gaussian prime factor, because each column has
both orientations among the rows. If an integer multiplicative
relation `product_i z_i^(c_i)=1` held, taking moduli would first give
`sum c_i=0`. Valuations at each `pi_S` then give `A^T c=0`.
Equation (1) forces `c=0`.

Hence actual primitive lattice configurations of arbitrarily many
points can be multiplicatively independent, not merely `B_k` for
one fixed `k`. This is **not** an endpoint counterexample: no claim
is made that these particular prime arguments place the points in
an arc of length `C sqrt(R)`.

## Implication for the remaining proof

The [all-combinations relaxation](integer_combination_separation.md)
uses precisely this full balanced incidence system, with compatible
short phases but without Gaussian-integral factors. Equation (1) makes
every nonzero zero-sum integer combination have nonzero combined
exponents. That construction additionally makes the corresponding
reduced monomial have imaginary part of absolute value at least one.
Its quotient by its conjugate is therefore not one. These two facts
together exclude every nontrivial equality of equal-length point
products, including lengths beyond thirty-two. Full exponent rank
alone would not suffice for arbitrary complex factor values, which
could have accidental multiplicative relations.

Thus neither a larger fixed product length nor a theorem excluding
all exact product relations would repair that relaxation. The
unresolved restriction must still couple the **actual Gaussian factor
values** to their **simultaneous short phases**. The finite width-four
result is useful, but full affine rank prevents extracting its
hypothesis merely by choosing a fixed number of these balanced rows.

This audit does not establish the sufficiency of the relaxed data for
lattice points. It identifies which arithmetic information the recent
product-rigidity results leave untested. The existing intrinsic
eight-point conductor inequality remains unproved.

Equation (1) and its row-subset consequence are elementary exact
linear algebra. The Gaussian realization uses distinct split primes;
its purpose is to verify that absence of relations alone is not a
cardinality bound. No Lean formalization is claimed.
