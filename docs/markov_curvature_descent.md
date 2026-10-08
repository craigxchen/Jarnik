# Directed-edge Markov curvature reduction

## Status

This note develops the graph/Markov framing using directed edges rather than
vertices.  It proves an entropy-versus-curvature inequality and shows that, at
the endpoint scale, all substantial branching entropy must be carried by the
evolution of common Gaussian divisors.  It does not yet prove the uniform arc
bound.

## 1. Directed-edge state space

Let `A` be a finite set of lattice points on the circle `|z|^2=N`, all lying in
an arc of Euclidean length at most `C N^(1/4)`.  A Markov state is a directed
edge `(i,j)` with `i != j`.  A transition chooses `k != i,j` and moves

    (i,j) -> (j,k).

Write

    v_ij = z_j-z_i.

The integer turning curvature is

    kappa(i,j,k) = |det(v_ij,v_jk)|.

Three distinct points on a circle are not collinear, so `kappa>=1`.

## 2. Exact multiplicity bound

Fix `(i,j)` and an integer `n>=1`.  The equation

    det(v_ij,z_k-z_j)=n

places `z_k` on one affine line parallel to `v_ij`; the equation with `-n`
places it on a second parallel line.  Each line intersects the circle in at
most two points.  Hence

    #{k : kappa(i,j,k)=n} <= 4.                       (2.1)

This is the basic non-flat Markov fact.

## 3. Entropy versus curvature

Let `K` be any random next vertex from the state `(i,j)`, and put

    X = kappa(i,j,K).

For every `s>1`, comparison with the probability distribution proportional to
`n^(-s)` gives

    H(K | i,j)
      <= log(4 zeta(s)) + s E[log X | i,j].           (3.1)

Equivalently,

    E[log X | i,j]
      >= (H(K|i,j)-log(4 zeta(s)))/s.                 (3.2)

Optimizing `s=1+1/max(H,1)` yields schematically

    E log X >= H(K|i,j)-log(H(K|i,j)+2)-O(1).         (3.3)

Thus high one-step branching entropy forces large integer curvature.

## 4. Endpoint geometric upper bound

If the three points lie in one endpoint arc, then the angle between the two
chord directions is at most the total angular width `O_C(N^(-1/4))`, while
each chord has length `O_C(N^(1/4))`.  Therefore

    kappa(i,j,k) <= C_1(C) N^(1/4),                  (4.1)

and

    log kappa <= (1/4) log N + O_C(1).               (4.2)

By itself this only gives a quarter-logarithmic entropy bound.

## 5. Removing the common Gaussian divisor

Let `gamma_ijk` be the maximal common Gaussian divisor of the three points,
and set

    G_ijk = Norm(gamma_ijk).

Both chord vectors are divisible by `gamma_ijk`, so

    G_ijk | kappa(i,j,k).

Define the reduced curvature

    rho(i,j,k) = kappa(i,j,k)/G_ijk in Z_{>=1}.       (5.1)

For a uniformly random triple from a nearly balanced squarefree cluster, each
conductor prime is common to all three with probability asymptotically at least
`1/4`.  Hence

    E log G_ijk >= (1/4) log N - o(log N).            (5.2)

Combining (4.2) and (5.2) gives

    E log rho(i,j,k) = o(log N)+O_C(1).               (5.3)

For an exactly balanced limiting process the logarithmic terms cancel
completely.

## 6. Conditional entropy decomposition

For fixed `(i,j)`, condition on the value of the triple common divisor `G`.
For each fixed pair `(G,n)`, the equation

    kappa = G n

again places the next vertex on at most two pairs of parallel lines, so there
are at most four possibilities.  Consequently, for any `s>1`,

    H(K | i,j)
      <= H(G_ijK | i,j)
         + log(4 zeta(s))
         + s E[log rho(i,j,K) | i,j].                (6.1)

This is the main Markov reduction.

At the endpoint, the last term is small on average.  Therefore any large
transition entropy must be carried by

    H(G_ijK | i,j),

namely by the choice of the new common-divisor state.

## 7. Interpretation as a killed descent process

The Markov chain should therefore be augmented by the current common Gaussian
divisor.  A transition either has large reduced curvature, which is expensive
by (6.1), or it moves into a new divisor state.  Dividing the corresponding
triple by that divisor reduces the circle norm and rescales all distances
exactly.

Thus the stochastic problem has been reduced to controlling a nested random
process of common-divisor states.  The missing theorem is an amortized divisor
entropy bound of the form

    sum_t H(G_{t+1} | current state) <= O_C(number of steps),

or an equivalent assertion that common-divisor choices cannot encode
unboundedly many independent bits while every reduced curvature remains small.

## 8. Exact obstruction

The remaining difficulty is not recurrence or cycle closure.  It is that a
pairwise common divisor may have many squarefree sub-divisors, and the next
vertex can choose among them.  The number of choices is not bounded by the
Euclidean geometry alone.

Any completion of the Markov proof must use compatibility between overlapping
states

    (i,j,k), (j,k,l), ...

to show that the common-divisor process loses conductor mass or entropy at a
uniform rate.  The directed-edge curvature theorem identifies this as the only
remaining source of large branching.