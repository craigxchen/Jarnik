# Weighted prime reveal, replica entropy, and recursive overlap

## Status

This note develops the prime-reveal/renormalization framework further. It proves
several exact identities, identifies one concrete obstruction to the naive
repair-support matching argument, and isolates a stronger recursive overlap
object that still retains the absolute cluster size.

## 1. Conductor-weighted random reveal

Assume the squarefree orientation model. Each split prime p contributes a
coordinate epsilon_p in {+1,-1} and arithmetic weight h_p = log p. Let
H = sum_p h_p = log N.

Instead of revealing coordinates in a fixed order, give each coordinate an
independent exponential clock of rate h_p. Equivalently, at each stage reveal
one remaining coordinate with probability proportional to h_p.

For two distinct successful leaves x,y, let D(x,y) be their set of differing
coordinates and d(x,y)=sum_{p in D(x,y)}h_p. The first revealed coordinate that
separates x and y then satisfies the exact formula

    P(first separator is p | x,y) = h_p / d(x,y)

for p in D(x,y). Endpoint separation gives d(x,y) >= H/2 - O_C(1), so this
filtration samples pairwise disagreement in the same metric used by the
Gaussian denominator bound.

## 2. Absolute entropy and replicas

Let A be the successful cluster, |A|=M, and sample X uniformly from A. At a
binary reveal-tree node, let the conditional branch masses be a and 1-a. Its
entropy contribution is

    h(a) = -a log a -(1-a)log(1-a).

There is an exact replica expansion

    h(a) = sum_{k>=1} [a(1-a)^k + (1-a)a^k]/k.

The kth summand is the probability that, after sampling X and k independent
comparison leaves Y_1,...,Y_k from the same conditional cluster, all comparison
leaves take the branch opposite to X, divided by k.

Summing over the probability-weighted reveal tree gives an exact representation
of log M by many-replica first-separation events. This is the scalar analogue
of preserving the successful branch weight in a live-coordinate martingale.

## 3. Arithmetic content of a replica event

Suppose coordinate p is the first revealed separator and X differs from all
Y_l at p. For every l, endpoint concentration gives a Gaussian near-relation

    2 s_p theta_p + sum_{q in R_l} 2 s_{l,q} theta_q = O_C(N^{-1/4}) mod 2pi,

where R_l uses only coordinates unrevealed at that time. The corresponding
conductor support has weight at least H/2-O_C(1).

Subtracting the l and m relations gives precisely the endpoint near-relation
between Y_l and Y_m. Hence a many-replica repair event is exactly a weighted
clique of successful leaves. It does not automatically produce independent
arithmetic relations beyond the pair relations already present in the cluster.

This is a concrete obstruction to the naive matching/core lemma: even if many
repair supports share the revealed prime p, the resulting relations may have
maximal algebraic dependence.

## 4. The stronger recursive object: completion profiles

At a reveal-tree node v, let f_v(t) be the number of successful completions
whose remaining angular contribution lies in the endpoint-scale cell t. If the
next prime has angle theta and children v+ and v-, then the exact recursion is

    f_v(t) = f_{v+}(t-theta) + f_{v-}(t+theta).

The cluster size at the root is the value of the final iterated operator at the
target interval. Balanced branching means both translated child profiles make
substantial contributions to the same target cell.

Thus entropy production is not merely the existence of a large repair support.
It is a quantitative overlap statement between two translated remaining-sum
profiles.

Define a scale-resolved overlap potential, schematically,

    Phi(f) = integral_0^infty || f/(f+u) ||_2^2 du,

or a finite discretized version. The scalar identity

    integral_0^infty (s/(s+u))^2 du = s

preserves total successful mass while spreading very small profile values over
all scales. A Bregman/resolvent inequality may then charge the change between
parent and child profiles to a probability-weighted entropy decrement without
an inverse-minimum-mass loss.

## 5. Renormalization target

Group primes into packets of comparable conductor height. Each packet acts on a
positive completion profile by convolution, translation along the successful
branch, and normalization. The target theorem is a contraction for the
scale-resolved overlap potential unless the profile is close to a coset
indicator of one common finite cyclic grid.

Exact localized fixed profiles are finite-grid coset indicators. Gaussian
unique factorization rules out a nontrivial exact persistent grid in the
primitive diffuse regime. What remains is quantitative stability of
approximate fixed profiles and, critically, persistence of the same grid over
many packet scales.

## 6. Current conclusion

The weighted reveal and replica identities are exact and align stochastic
separation with conductor weight. They also show why the original
repair-support hypergraph is insufficient: replica repairs can be completely
dependent only in appearance, because all their differences are already the
existing clique relations.

The correct renormalized state is therefore the full positive completion
profile, not a selected repair support. The next load-bearing estimate is a
scalar resolvent-Bregman inequality for the two-child recursion, followed by an
inverse theorem for near-zero dissipation under repeated Gaussian packet
convolutions.
