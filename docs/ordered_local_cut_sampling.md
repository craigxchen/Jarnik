# Ordered local cut sampling from nonpositive pair correlations

This is a quantitative ordered-average application of the obtuse Bessel
inequality in [uniform profile extraction](uniform_profile_extraction.md),
not a new extraction theorem. That note already extracts one tuple with
all pattern weights almost uniform; sorting such a tuple preserves that
property. Here the conclusion controls the average absolute discrepancy
over **all uniformly sampled subsets in their prescribed order**. It
applies to order-sensitive local tests, without exchangeability of the
rows or layers. It does not apply to consecutive windows or to an
arbitrarily biased selection of subsets.

## Setup and ordered-average theorem

Let `s_1,...,s_M` be sign functions on a probability space of layers, with

```text
E_layers(s_i s_j) <= 0  for i != j.                    (1)
```

These are raw correlations, not centered covariances. A finite layer
system of positive total weight `W` is normalized by dividing its weights
by `W`. Fix any order on the rows. For a uniform `k`-subset
`I={i_1<...<i_k}`, put

```text
P_I(f) = E_layers f(s_i1,...,s_ik),
fhat_A = 2^(-k) sum_(epsilon in {+1,-1}^k) f(epsilon) epsilon_A,
fbar = fhat_empty,
||fhat_r||_2^2 = sum_(|A|=r) fhat_A^2.
```

For every real function `f` on the `k`-cube,

```text
E_I |P_I(f)-fbar|
 <= sum_(r=1..k) ||fhat_r||_2 sqrt(2 binom(k,r)/(M-r+1)).   (2)
```

For fixed `k,f`, the right side is `O_f(M^(-1/2))`. If `f` is
complement-invariant, its odd Fourier degrees vanish. Formula (2) is
valid for every prescribed permutation of the rows, including their
geometric order. No permutation is averaged in its statement.

For the indicator of any specified ordered sign pattern `epsilon`, and
for the indicator of that pattern together with its complement,
respectively, (2) implies

```text
E_I |P_I(epsilon)-2^(-k)|
 <= (1-2^(-k)) sqrt(2/(M-k+1)),

E_I |P_I({epsilon,-epsilon})-2^(1-k)|
 <= (1-2^(1-k)) sqrt(2/(M-k+1)).                         (3)
```

Thus most subsets, not just one extracted subset, have small discrepancy
for any fixed test: Markov's inequality bounds the exceptional fraction
by the right side of (2) divided by the desired tolerance.

## Counting proof: rank positions must be summed

For a row set `J`, write `mu_J=E_layers product_(j in J) s_j`.
The existing obtuse Bessel inequality is

```text
sum_i <s_i,h>^2 <= 2 ||h||_2^2.
```

Apply it with `h` the product on an `(r-1)`-set, discard terms whose
index lies in that set, and sum over all such sets. Each `r`-set occurs
exactly `r` times, so

```text
sum_(|J|=r) mu_J^2
 <= (2/r) binom(M,r-1)
 = [2/(M-r+1)] binom(M,r).                            (4)
```

A fixed subset of **rank positions** in a sorted sample does not, in
general, sample row sets uniformly. The correct identity sums over all
rank-position subsets of a fixed size:

```text
E_I sum_(|A|=r) mu_(I_A)^2
 = [binom(M-r,k-r)/binom(M,k)] sum_(|J|=r) mu_J^2
 <= 2 binom(k,r)/(M-r+1).                             (5)
```

Indeed each `J` is contained in `binom(M-r,k-r)` samples `I`, and in each
sample it occupies exactly one rank-position subset `A`. Fourier
expansion, the triangle inequality, Cauchy--Schwarz in `A`, and then
Cauchy--Schwarz in `I` give (2). For a pattern indicator all degree-`r`
coefficients have magnitude `2^(-k)`; for a pattern/complement pair all
even-degree coefficients have magnitude `2^(1-k)` and odd ones vanish.
Summing their binomial counts and replacing each denominator by `M-k+1`
gives (3).

## A direct block bound for complement-invariant tests

There is also an elementary bound useful when the Fourier `l^1` norm is
small. Partition the row order into `K>=k` consecutive nonempty blocks,
of sizes between `n_min` and `n_max`. Set

```text
L_f = sum_(A nonempty) |fhat_A|,
osc(f) = max f - min f.
```

For complement-invariant `f`,

```text
|E_I P_I(f)-fbar|
 <= L_f/n_min
    + osc(f) binom(k,2)(n_max-1)/(M-1).                (6)
```

For block `b`, its sign bias `x_b=n_b^(-1) sum_(i in b) s_i` satisfies
`E_layers x_b^2<=1/n_b` by (1). Conditional on the sample meeting
distinct blocks, and on those blocks, its selected rows are independent
and uniform within their respective blocks. The conditional Fourier
formula therefore replaces each sign by the corresponding `x_b`.
Every nonconstant monomial of a complement-invariant test has at least
two factors. Bounding all but two in absolute value by one and applying
Cauchy--Schwarz bounds its layer expectation by `1/n_min`. Finally, the
probability of a block collision is at most
`binom(k,2)(n_max-1)/(M-1)`, by a union bound on sampled pairs. On that
event the discrepancy from `fbar` is at most `osc(f)`. This proves (6).

With `K=floor(sqrt(M))>=max(k,2)` and almost equal block sizes,
`n_min>=K` and `(n_max-1)/(M-1)<=1/(K-1)`. This gives explicitly
`L_f/K + osc(f) binom(k,2)/(K-1)`. Unlike (2), (6) is only a bound on
the signed average discrepancy.

## Arithmetic scope and the ordered five-point test

For the test `01110` or `10001`, `fbar=1/16` and `L_f=15/16`.
Consequently its average inherited layer weight is `(1/16+o(1))W` in
any ordered sign system satisfying (1). The same statement holds for
the two constant patterns. The
[ordered conductor-index criterion](ordered_conductor_index_uniformity_criterion.md)
uses sharper special estimates and the actual integer index divisibility
theorem. Removing a selected tuple's common Gaussian divisor removes
exactly its constant threshold layers. Thus the two limiting fractions
are `1/16` for the distinguished cut and `15/16` for its intrinsic total
conductor. This explains the sufficient limiting intrinsic exponent
`delta<1/15` when that arithmetic criterion charges `delta log N_I`.

The present theorem supplies sampling only. It proves neither the
required arithmetic upper bound on the conductor-supported index nor
a uniform circle-point bound. In particular, orientation invariance of
a cut test does not itself supply any Gaussian-integrality inequality.

The accompanying [exact checker](check_ordered_local_cut_sampling.py)
tests ordered Fourier counting, degreewise discrepancy bounds, and the
block estimate on permuted and individually column-oriented Hadamard
fixtures. The proof above, rather than those finite checks, establishes
the general result.
