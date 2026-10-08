# Markov triple-gcd entropy correction

## Status

This note audits the proposed directed-edge Markov proof. It proves an exact
state-conditioned factorial lower bound for reduced curvatures and shows that
the previously proposed low-entropy theorem for triple common-divisor states is
false in the required strength. The uniform endpoint theorem is not proved.

## Setup

Let `z_1,...,z_M` be lattice points on one circle of norm `N`, contained in an
endpoint arc. Fix a directed edge `(i,j)` and write

    kappa_k = |det(z_j-z_i, z_k-z_j)|

for `k != i,j`. The geometry of the endpoint arc gives

    kappa_k <= C^3 N^(1/4).

Let `gamma_k` be a maximal Gaussian divisor common to `z_i,z_j,z_k`, and put

    G_k = Norm(gamma_k),
    rho_k = kappa_k/G_k.

Then `rho_k` is a positive integer.

## Bounded multiplicity at fixed state

For fixed values `G` and `rho`, the condition

    kappa_k = G rho

puts `z_k` on one of two parallel lattice lines determined by the fixed chord
`z_j-z_i`. Each line meets the circle in at most two points. Hence

    #{k : (G_k,rho_k)=(G,rho)} <= 4.

If `n_G=#{k:G_k=G}`, it follows that

    product_{k:G_k=G} rho_k >= (floor(n_G/4)!)^4.

Summing logarithms and using Stirling gives

    sum_k log rho_k
      >= (M-2) log(M-2)
         -(M-2) H(p_G)
         -O(M),

where `p_G=n_G/(M-2)` and `H(p_G)` is the Shannon entropy of the triple-gcd
state distribution for the uniform continuation from `(i,j)`.

Equivalently,

    H(p_G)
      >= log(M-2)
         -(1/(M-2)) sum_k log rho_k
         -O(1).

Thus small reduced curvature does not force low state entropy. It forces the
opposite: the triple-gcd state must nearly identify the next vertex.

## Arithmetic upper bound

Since `rho_k=kappa_k/G_k` and `kappa_k <= C^3 N^(1/4)`, one has

    sum_k log rho_k
      <= (M-2)(1/4 log N + O_C(1))
         -sum_k log G_k.

For a nearly balanced squarefree cluster, the primewise calculation of
`sum_k log G_k` is near `(M-2) log N/4` for a typical nearly critical edge.
Consequently the upper bound for `sum log rho_k` can be only `o(M log N)+O_C(M)`.
Inserted into the factorial inequality, this says that a large cluster forces

    H(p_G) = log M - o(log M)

whenever the reduced-curvature budget is small.

## Failure of the proposed finish

The hoped-for theorem was that overlapping directed edges make the process
`G_k` low entropy. This is false at the combinatorial level. For a fixed edge,
`G_k` records the subset of primes on which the third vertex agrees with the two
edge endpoints. In balanced binary orientation systems this subset can encode
`k` almost injectively. Consecutive triple-gcd states are not nested, so edge
overlap does not produce an amortized `O(1)` entropy bound.

Therefore the directed-edge Markov framework currently proves the following
rigidity statement, not a uniform bound:

> If many endpoint points coexist and reduced integer curvatures remain small,
> then for each nearly critical edge the triple common-divisor state must encode
> the continuation vertex with bounded ambiguity.

Any successful continuation must exploit the arithmetic realizability of this
near-injective divisor encoding. Entropy or Markov overlap alone is insufficient.
