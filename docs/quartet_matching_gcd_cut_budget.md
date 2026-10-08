# Quartet matching content: a global gcd theorem and the fair-cut obstruction

For four distinct Gaussian circle points, let `gamma_Q` be the Gaussian
gcd of the three products of differences associated with their perfect
matchings, and put `Gamma_Q=Norm(gamma_Q)`. Any two of those products
generate the same ideal, because the third is their difference, up to
sign. Thus this definition is independent of ordering and agrees with
the factor `Pi/s^2` in
`moving_quadrilateral_fifth_point_product.md`.

This note proves an exact global-gcd description, a nested prime-cut
divisibility statement, and a fair-profile obstruction. The obstruction
can be realized by actual primitive equal-norm Gaussian rows in arcs of
angular width tending to zero. No endpoint bound on those arcs is claimed.
In particular, neither primitivity nor a small unnormalized angle makes
some quartet content bounded. A consequence at the endpoint scale would
have to use more than the valuation budgets calculated here.

## 1. Exact global gcd of all quartet contents

For arbitrary distinct Gaussian integers `z_1,...,z_m`, `m>=4`, define

\[
\mathfrak g=\gcd_{|Q|=4}\gamma_Q.
\]

Equivalently, this is the gcd of all products
`(z_i-z_j)(z_k-z_l)` with four distinct indices. Fix a Gaussian prime
`pi`, put

\[
r=\min_{i<j}v_\pi(z_i-z_j),
\]

and partition the indices by the residues of `(z_i-z_1)/pi^r` modulo
`pi`. There are at least two classes. If no class has `m-1` indices,
then

\[
\boxed{v_\pi(\mathfrak g)=2r.}\tag{1}
\]

If a class `B` has `m-1` indices, put
`r_B=min_{i<j in B} v_pi(z_i-z_j)>r`. Then

\[
\boxed{v_\pi(\mathfrak g)=r+r_B.}\tag{2}
\]

For (1), the complete multipartite graph of pairs from different classes
has a matching of size two: its largest class has size at most `m-2`
and `m>=4`. The two edge valuations are `r`, attaining the universal
lower bound `2r`. For (2), every four-index matching either has one edge
from the sole outlier to `B` and one edge within `B`, or has both edges
within `B`. Its valuation is at least `r+r_B`. A pair attaining `r_B`,
and an outlier edge using a third member of `B`, attains this bound.

Thus global content is governed by a very specific near-total residue
concentration, after the common difference factor is removed.

For points on one circle there is a further conjugation identity:

```text
Norm(gcd_Q gamma_Q)=gcd_Q Norm(gamma_Q).
```

The [independent audit, Section 6](quartet_matching_gcd_cut_budget_independent_audit.md)
proves simultaneous attainment of the two conjugate valuation minima at
each rational prime. In a primitive tuple, the generic witnessing quartet
uses the two smallest and two largest valuation levels; an all-but-one
extreme case uses a minimum-residue majority pair. The quartet can vary
with the prime. Thus the identity does not turn the common gcd into the
minimum of the individual norms.

There is no primitive-circle theorem saying `g` divides two. For any
`e>=1`, the five points

\[
\pi^e,\quad i\pi^e,\quad-\pi^e,\quad-i\pi^e,\quad\bar\pi^e,
\qquad \pi=2+i,
\]

are distinct, have common norm `5^e`, and have Gaussian gcd one.
The first four form the unique four-element zero residue class modulo
`pi`. Their within-class minimum difference valuation is `e`, so (2)
gives `v_pi(g)=e`. In particular the odd part of the global content can
be arbitrarily large. These five points are not a short-arc example.

## 2. Nested cut contribution and an all-quartet divisibility

Now let all `m` points have norm `N=R^2` and Gaussian gcd one. Then `N`
is odd and all its prime factors split. Fix `p^e || N`, choose `pi`
above `p`, and put `t_i=v_pi(z_i)`. For a quartet write its four values
in increasing order `t_(1)<=...<=t_(4)`. Matching minima give

\[
v_p(\Gamma_Q)\ge
2e+t_{(1)}+t_{(2)}-t_{(3)}-t_{(4)}.
\tag{3}
\]

Indeed a matching's two smaller endpoint valuations have sum at least
the two smallest `t` values. Pairing small with large attains that
combinatorial minimum; equal-level cancellations can only raise the
actual difference valuations. Apply the same statement at the conjugate
prime, where the levels are `e-t_i`, to obtain (3).

For `1<=h<=e`, put `S_h={i:t_i>=h}`. Equation (3) is equivalently

\[
\boxed{v_p(\Gamma_Q)\ge
\sum_{h=1}^e\bigl||Q\cap S_h|-2\bigr|.}\tag{4}
\]

This handles arbitrary nested levels. It does not assume that primes
occur only to exponent one or that the cuts are independent. The exact
formula with the additional residues is equation (16) of the preceding
moving-quadrilateral note.

For a cut of size `s`, summing the layer contribution over all quartets
gives exactly

\[
J_m(s)=2\left[\binom s4+\binom{m-s}4\right]
 +(m-s)\binom s3+s\binom{m-s}3.\tag{5}
\]

Its minimum is at the balanced cut. To check this directly, put
`t=m-s-1`. Moving one label to the `s` side gives

\[
J_m(s+1)-J_m(s)
=\frac{(s-t)((m-2)(m-3)+2st)}6.\tag{6}
\]

The sign is therefore the sign of `2s+1-m`. Write
`J_m=J_m(floor(m/2))`. Summing (4), including all prime factors of `N`,
proves the exact integer divisibility

\[
\boxed{N^{J_m}\mid\prod_{|Q|=4}\Gamma_Q.}\tag{7}
\]

This is a lower product bound. It does not force any individual
`Gamma_Q` to be small. For comparison, on an endpoint arc of length
`C sqrt(R)`, the immediate matching bound is
`Gamma_Q<=C^4 N`, whereas (7) has mean exponent approaching `3/4`.
These two estimates do not contradict one another.

## 3. Exact comparison with triangle content

For every quartet, the preceding note proves

\[
\Gamma_Q=\frac{\prod_{I\subset Q,\ |I|=3}|D_I|}{s_Q^2},
\qquad s_Q\in\mathbb Z_{>0}.
\]

Thus globally

\[
\boxed{\prod_Q\Gamma_Q
=\frac{(\prod_{|I|=3}|D_I|)^{m-3}}{\prod_Qs_Q^2}.}\tag{8}
\]

At a simple two-level cut with no excess same-level residues, a triangle
has determinant valuation `e` when monochromatic and zero otherwise.
A quartet's `s_Q` has valuation `e` when monochromatic and zero otherwise;
its `Gamma_Q` has valuation `e|k-2|`, where `k` is its cut intersection.
Consequently the exact accounting identity is

\[
(m-3)\left[\binom s3+\binom{m-s}3\right]
=J_m(s)+2\left[\binom s4+\binom{m-s}4\right].\tag{9}
\]

The quartet product has redistributed the triangle product's pure-cut
weight into matching content and `s_Q`; it has not created an extra
conductor exponent. In this residue-free model the triangle gcd has
zero valuation whenever the cut is nontrivial, so this comparison also
applies to the primitive triangle data of
`affine_shape_radius_divisibility.md`.

## 4. Every quartet can have large content under fair balanced cuts

Let `m=2n>=6`, and give every `n`-element cut equal logarithmic norm
weight `w`. There are `K=binom(m,n)` such cuts. Each fixed quartet has
the same sum of weights in (4), by permutation symmetry. Its normalized
mean contribution is

\[
\boxed{\beta_m=\frac{J_m}{\binom m4}
=\frac{3(n-1)(n-2)}{(2n-1)(2n-3)}
\longrightarrow\frac34.}\tag{10}
\]

For example `beta_6=2/5`, `beta_8=18/35`, and `beta_10=4/7`.
Hence a profile with these weights forces, simultaneously for every
quartet,

\[
\boxed{\log\Gamma_Q\ge(\beta_m-o(1))\log N.}\tag{11}
\]

Nevertheless every conductor cut has at least two labels on each side.
There is a quartet meeting it in two labels on each side; pairing across
the cut gives valuation zero at both conjugate primes, without any
cancellation issue. Therefore the global gcd of all `gamma_Q` has
**zero valuation at every conductor prime** in this profile. Coprimality
of their common conductor support does not make one quartet content
bounded: different quartets retain different parts of the conductor.
Extra primes from accidental congruences outside `N` are not being
declared absent.

## 5. Actual primitive rows realizing the fair profile

The preceding lower-bound profile is compatible with actual Gaussian
points, primitivity, absence of collinear triples, and small angular
width. Here is an elementary construction, without any assumption about
Gaussian prime angles.

Enumerate the `K` balanced cuts by `j=1,...,K`. Choose a positive integer
`M` divisible by every nonzero `j^2-k^2`, and set, for integers `q>=1`,

\[
\kappa_j=2jMq+i,\qquad n_j=\operatorname{Norm}\kappa_j.
\]

The odd integers `n_j` are pairwise coprime. If a prime divided both
`4j^2M^2q^2+1` and `4k^2M^2q^2+1`, it would not divide `2Mq`, but
their difference would force it to divide `j^2-k^2`, and hence `M`, a
contradiction. Also `gcd(kappa_j,conj(kappa_j))=1`, since a common
Gaussian prime would divide `2i` while `n_j` is odd.

For the cut `S_j`, define the row

\[
z_i=\prod_{j:i\in S_j}\kappa_j
     \prod_{j:i\notin S_j}\bar\kappa_j.
\tag{12}
\]

These rows all have norm `N=product_j n_j`. Their Gaussian gcd is one,
because every block is allocated to both orientations and distinct
blocks have disjoint prime support. They are distinct: some balanced cut
separates any chosen pair of labels, and its nonunit prime factors have
different allocations in those rows. No three distinct points on a
circle are collinear.

Every prime-power factor of a block has exactly its prescribed balanced
two-level allocation, with the conjugate orientation possibly reversed.
Moreover

\[
\log n_j=2\log q+O_m(1),\qquad
\log N=2K\log q+O_m(1).
\]

Thus (11) holds for every quartet of these actual configurations.
The angular width is in fact comparable to `1/q`, not just bounded above
by it. Write `sigma_ij=1` or `-1` for the chosen orientation and put
`h_i=sum_j sigma_ij/j`. The unique largest power of two among
`1,...,K` supplies the unique summand of smallest 2-adic valuation in
each `h_i`; consequently every `h_i` is nonzero. Balanced columns give
`sum_i h_i=0`, so these values are not all equal. The lifted arguments
have the expansions

```text
arg(z_i)=h_i/(2Mq)+O_m(q^(-3)).
```

Their span is therefore comparable to `1/q`, with positive constants
depending on m. Since `R` has order `q^K`, the endpoint-normalized span
is comparable to `q^(K/2-1)` and tends to infinity. The construction is therefore
an obstruction to the valuation-level selection argument, including
its primitive-row and unscaled angular premises. It is not a family of
endpoint counterexamples and does not rule out a selection theorem
that uses the precise endpoint phase accuracy.

The identities and local matching formula are checked in
`check_quartet_matching_gcd_cut_budget.py`. They provide a concrete
quartet theorem and a precise limit on extracting a bounded quartet
content from fair cut averages alone.
