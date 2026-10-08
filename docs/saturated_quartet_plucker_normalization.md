# Exact deletion factors and saturated integral Plücker coordinates

Let `z_0,...,z_(m-1)` be distinct Gaussian integers, with `m>=4`. They
need not lie on a circle for the algebraic statements in Sections 1–3.
Let

\[
D=\gcd_{i<j}(z_j-z_i),\qquad
u_i=(z_i-z_0)/D\in\mathbb Z[i].
\]

Thus `u_0=0` and the Gaussian gcd of all `u_i-u_j` is a unit. Define
the deletion factors

\[
A_i=\gcd_{j<k,\ j,k\ne i}(u_k-u_j),\qquad A=\prod_i A_i.
\tag{1}
\]

If `G` is the Gaussian gcd of all products of two disjoint original
differences, then

\[
\boxed{\gcd(A_i,A_j)=1\ (i\ne j),\qquad G\sim D^2A.}\tag{2}
\]

Here `~` means equality up to a Gaussian unit. Moreover

\[
\boxed{e_{ij}=\frac{A_iA_j}{A}(u_j-u_i)\in\mathbb Z[i]}\tag{3}
\]

are primitive integral Plücker coordinates: their common Gaussian gcd
is one and they are the two-by-two minors of an explicit integral
two-by-`m` matrix. Each column of that matrix is itself primitive.
This latter assertion uses the exact local deletion factors, not merely
the fact that the gcd of all minors is one.

## 1. The deletion factors are exactly the remaining matching content

Fix a Gaussian prime `pi`. Because the normalized differences have gcd
one, at least two residues occur among the `u_i` modulo `pi`.

A prime can divide `A_i` only if all rows other than `i` have the same
residue. The outlier `i` must have a different residue, since otherwise
the prime would divide every normalized difference. Therefore there
is at most one index `i` for which `pi|A_i`, proving pairwise
coprimality.

If there is no such index, no residue class has `m-1` members. Two
disjoint pairs from different classes then give a matching product
with valuation zero. If there is a unique outlier `i`, write

\[
v_\pi(A_i)=a
=\min_{j<k,\ j,k\ne i}v_\pi(u_k-u_j)>0.
\]

Every two-edge matching either uses an outlier edge and an edge within
the majority, or two within-majority edges. Its valuation is at least
`a`. A minimizing majority pair, together with an edge from the outlier
to a third majority member, attains `a`; that third member exists
because `m>=4`. Hence the gcd of normalized matching products has
valuation exactly `a`, proving (2) prime by prime.

Equivalently, for each edge `(i,j)`, every `A_k` with `k` different
from both endpoints divides `u_j-u_i`. Their pairwise coprimality
gives the integrality in (3).

## 2. Primitive minors and primitive columns

Extend `e_ij` antisymmetrically and put `e_ii=0`. They obey both

\[
e_{ij}e_{kl}-e_{ik}e_{jl}+e_{il}e_{jk}=0,
\qquad
A_i e_{jk}-A_j e_{ik}+A_k e_{ij}=0.\tag{4}
\]

These follow directly from (3) and the corresponding identities for
differences. To check primitivity without assuming it from the
Plücker equations, fix a Gaussian prime.

* If the prime divides no `A_i`, formula (3) changes every edge
  valuation by zero. Some original normalized edge has valuation zero.
* If it divides `A_i` to exponent `a`, the prime divides no other
  `A_j`. The exceptional-residue description in Section 1 gives
  `v_pi(u_j-u_i)=0` for every `j!=i`. Formula (3) therefore gives
  `v_pi(e_ij)=a+0-a=0` for every such `j`.

Thus the common gcd of all `e_ij` is one.

Choose Gaussian Bezout coefficients with `sum_i c_i A_i=1`, and define

\[
b_j=\sum_i c_i e_{ij}\in\mathbb Z[i].
\tag{5}
\]

The weighted three-index identity in (4) gives

\[
\boxed{e_{ij}=A_i b_j-A_j b_i.}\tag{6}
\]

Indeed `A_i e_kj-A_j e_ki=A_k e_ij`; multiply by `c_k` and sum.
This constructs the promised matrix with columns `(A_i,b_i)^T`.
Because its minors have gcd one, its image is all of `Z[i]^2`, or
equivalently it is a saturated rank-two realization.

For a prime not dividing `A_i`, the `i`th column is primitive at that
prime. If `pi|A_i`, Section 2's second case gives `v_pi(e_ij)=0` for
every `j!=i`, and `A_j` is a unit modulo `pi`. Equation (6) then
forces `b_i` to be a unit modulo `pi`. Consequently

\[
\boxed{\gcd(A_i,b_i)=1\quad\text{for every }i.}\tag{7}
\]

This proof addresses the individual columns explicitly: a primitive
matrix of minors alone would not establish (7).

## 3. What the normalization retains

Put `U=sum_i c_i A_i u_i`, a Gaussian integer. Equations (3) and (5)
give the exact reconstruction

\[
\boxed{b_j=\frac{A_j}{A}(u_j-U),\qquad
u_j=U+A\frac{b_j}{A_j}.}\tag{8}
\]

Thus the primitive columns describe rational points `b_j/A_j`, with
pairwise coprime reduced denominators. Both the denominators and the
common scale `A` are retained. If the original rows lie on an
origin-centered circle, these rational points lie on a translated
circle of radius `R/(|D A|)` and center
`-(z_0/D+U)/A`; equation (8), rather than an integral-center assertion,
is the normalization statement.

In the primitive endpoint setting of the preceding content theorem,

\[
\operatorname{Norm}A
=\frac{\operatorname{Norm}G}{(\operatorname{Norm}D)^2}
\le\frac{C^6}{32(\operatorname{Norm}D)^2}R^{2\gamma_m}
\quad(m\ge9).
\]

This is a restriction on the total denominator product. It does not
bound the integral numerators `b_i`, and it does not select a quartet
of bounded individual matching content. No cardinality conclusion
beyond the earlier content bound is inferred.

## 4. Exact constant-content circle configurations of growing size

Here is a concrete obstruction to using saturated Plücker integrality,
circle identities, and a small *unnormalized* angle as sufficient
replacement hypotheses for endpoint geometry.

Choose `k>=2` Gaussian blocks

\[
\kappa_j=x_j+i,\qquad x_j>0\text{ even},
\]

whose norms `n_j=x_j^2+1` are pairwise coprime. For every sign vector
`epsilon in {+1,-1}^k`, take

\[
z_\epsilon=\prod_{j=1}^k(x_j+\epsilon_j i).
\tag{9}
\]

There are `m=2^k` distinct points, all of norm `N=product_j n_j`,
and their Gaussian gcd is one. Each block and its conjugate are
coprime, and distinct blocks have disjoint support, proving both
primitivity and distinctness by their allocations.

Every difference is divisible by `2i`. Conversely, differences along
an edge of the sign cube, obtained by changing only sign `j`, are

\[
2i\prod_{l\ne j}(x_l+\epsilon_l i),
\]

up to sign. The gcd of the products on the right, over all choices
of the other signs, is one. Hence exactly `D~2i`.

Choose two disjoint parallel edges changing sign `j`, with opposite
choices at every other coordinate. Their difference product is

\[
-4\prod_{l\ne j}n_l=-4N/n_j
\]

up to sign. The integer gcd of these quantities as `j` varies is
four, because the `n_j` are pairwise coprime. Meanwhile every matching
product is divisible by `(2i)^2=-4`. Therefore exactly

\[
\boxed{G\sim4,\qquad \operatorname{Norm}G=16,\qquad
A_i\text{ is a unit for every }i.}\tag{10}
\]

This gives saturated integral Plücker configurations of arbitrarily
large cardinality with absolutely bounded global matching content,
on actual primitive origin circles. There are no collinear triples.

One explicit choice with arbitrarily small angular width is
`x_j=2jMq`, where `M` is divisible by all nonzero `j^2-l^2` and
`q` tends to infinity. The norms are pairwise coprime: a common prime
would divide `j^2-l^2`, hence `M`, contradicting that it divides
`4j^2M^2q^2+1`. The angular span is exactly
`2 sum_j arctan(1/x_j)` once that sum is less than `pi/2`, and is of
order `1/q`. The radius has order `q^k`. For `k>=3`, the endpoint
constant therefore grows as `q^(k/2-1)`.

Thus this is not an endpoint counterexample. It shows precisely that
constant global content and the full integral Grassmannian realization
do not themselves force a cardinality bound; the missing information
is the endpoint angular accuracy, not a missing Plücker identity.

The checker `check_saturated_quartet_plucker_normalization.py` verifies
the factorization, integral primitive minors, explicit Bezout
realization, primitive individual columns, reconstruction, and the
constant-content sign-cube examples with exact Gaussian arithmetic.
It passes on 600 arbitrary tuples, 100 tuples with deliberately forced
deletion factors, and four primitive 8- and 16-point circle examples.
