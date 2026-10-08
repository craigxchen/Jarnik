# Fixed-edge divisor-state normalization

## Status

This note continues the directed-edge Markov formulation. It proves that, in a
nearly balanced squarefree clique, the triple common-divisor state attached to a
continuation of one fixed edge has norm at the quarter-conductor scale, and that
conditioning on this state reduces the three points to the classical lattice-arc
scale `R^(1/3)`. Hence each state has uniformly bounded multiplicity.

It does not bound the number of distinct states. That remaining counting problem
is equivalent to controlling near-middle Gaussian divisors of one fixed pairwise
common factor, with additional critical-arc realizability constraints.

## 1. Triple unanimous weight

Let the squarefree conductor be `N`, with `H = log N`. Represent each point by
its Gaussian orientation sign vector. For three points `i,j,k`, let `D_ab` be
the conductor-weighted Hamming distance between points `a,b`, and let `U_ijk`
be the total conductor weight of primes at which all three orientations agree.

At each prime, either all three signs agree, contributing zero to

    D_ij + D_jk + D_ki,

or one sign is exceptional, contributing twice the prime weight. Therefore

    D_ij + D_jk + D_ki = 2(H - U_ijk),

and hence

    U_ijk = H - (D_ij + D_jk + D_ki)/2.              (1.1)

For a fixed balanced clique extracted from a hypothetical unbounded family,
all pair distances satisfy

    D_ab = H/2 + o(H).

Thus

    U_ijk = H/4 + o(H).                              (1.2)

If `gamma_ijk` is the maximal oriented Gaussian divisor common to all three
points, then

    log Norm(gamma_ijk) = U_ijk,

so

    Norm(gamma_ijk) = N^(1/4+o(1)).                  (1.3)

## 2. Division reaches the classical one-third scale

The original circle has radius

    R = N^(1/2)

and the endpoint arc has length

    L <= C N^(1/4).

Divide the three points by `gamma = gamma_ijk`. The new circle has squared norm

    N' = N / Norm(gamma) = N^(3/4+o(1)),

so its radius is

    R' = N^(3/8+o(1)).

Euclidean lengths are divided by

    |gamma| = Norm(gamma)^(1/2) = N^(1/8+o(1)).

Thus the descended arc length is

    L' <= C N^(1/8+o(1)) = C (R')^(1/3+o(1)).        (2.1)

At exact balance, the exponent is precisely `1/3`. Classical lattice-arc
geometry therefore gives a bound depending only on `C` (and on a fixed balance
margin) for the number of descended lattice points associated with one fixed
common-divisor state.

This recovers the bounded-multiplicity conclusion from the curvature argument
in a more geometric form.

## 3. Fixed-edge state space

Fix a nearly critical edge `(i,j)`. Let `g_ij` be the maximal Gaussian divisor
common to its endpoints. Since

    D_ij = H/2 + o(H),

we have

    Norm(g_ij) = N^(1/2+o(1)).                       (3.1)

For each continuation `k`, the triple state `gamma_ijk` divides `g_ij`, and by
(1.3)

    Norm(gamma_ijk) = N^(1/4+o(1)).                  (3.2)

Thus every continuation chooses a near-middle Gaussian divisor of the fixed
Gaussian integer `g_ij`.

The quotient

    delta_ijk = g_ij / gamma_ijk

also has norm `N^(1/4+o(1))`. Hence the state is equivalently a near-balanced
factorization

    g_ij = gamma_ijk delta_ijk                       (3.3)

into two Gaussian factors of comparable norm.

Conditioning on `gamma_ijk` leaves only uniformly many possible continuations,
but the number of possible factorizations (3.3) can grow with the number of
split prime factors of `g_ij`.

## 4. Exact Markov interpretation

For the directed-edge chain, the continuation entropy decomposes as

    H(K | i,j)
      = H(gamma_ijK | i,j)
        + H(K | gamma_ijK, i,j).                     (4.1)

The second term is `O_C(1)` by the one-third-scale reduction. Therefore

    H(K | i,j) = H(gamma_ijK | i,j) + O_C(1).        (4.2)

So the Markov process does not merely suggest that the divisor state is
important: asymptotically, all branching entropy is exactly the entropy of the
near-middle factorization of `g_ij`.

## 5. Remaining finite arithmetic problem

A uniform endpoint bound would follow from a uniform bound on the number of
near-middle divisors `gamma | g_ij` that are realizable as maximal triple gcds
with a third point in the endpoint arc.

The unrestricted divisor count is not uniformly bounded. A squarefree Gaussian
integer with many prime factors can have many divisors of norm close to the
square root of its norm. Therefore the required theorem must use the additional
realizability condition:

> both endpoints divided by `gamma`, together with a primitive third point,
> lie on one circle in an arc of critical one-third-scale length, and `gamma`
> is the maximal common Gaussian divisor of the triple.

This is a concrete fixed-edge divisor-realizability problem. It is the exact
load-bearing statement left by the Markov formulation.

## 6. Warning

The overlap of consecutive directed edges does not make the states
`gamma_ijk` nested. Abstract balanced sign systems can make the map

    k -> gamma_ijk

nearly injective for every fixed edge. Thus no purely entropy-theoretic or
Boolean argument can bound the state count. Any completion must exploit the
Gaussian arguments or the critical one-third-scale geometry after division.
