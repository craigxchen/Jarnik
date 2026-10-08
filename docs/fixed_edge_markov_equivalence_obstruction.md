# Fixed-edge Markov divisor states: exact obstruction

## Status

This note closes the current Markov/divisor-state route. It does not prove the
uniform endpoint bound. It shows that the proposed fixed-edge
"divisor-realizability theorem" is not a genuine simplification: after the
state normalization, different divisor states descend to unrelated circles,
and a uniform bound on the number of realizable states is essentially
Equivalent to the original endpoint theorem.

## 1. Fixed edge and state decomposition

Let `z_i,z_j` be two points of norm `N`, and let `g` be their maximal common
Gaussian conductor divisor. Write

    z_i = g a,
    z_j = g b,

so

    Norm(a)=Norm(b)=N/Norm(g).

For a third point `z_k`, let

    gamma = gcd_Z[i](z_i,z_j,z_k)

and write

    g = gamma delta,
    z_k = gamma c.

Then

    z_i/gamma = delta a,
    z_j/gamma = delta b,
    z_k/gamma = c

all lie on the circle

    Norm(w)=N/Norm(gamma).

The original arc is scaled by `1/|gamma|`.

The maximality of `gamma` says that `c` has no further Gaussian conductor
factor common with both `delta a` and `delta b`.

## 2. Different states do not share a descended problem

The descended radius is

    R_gamma = sqrt(N/Norm(gamma)),

and the descended arc length is

    L_gamma = L/sqrt(Norm(gamma)).

Both depend on `gamma`. Moreover the fixed descended endpoints are

    delta a, delta b,

and `delta=g/gamma` also depends on `gamma`.

Thus two different states `gamma_1,gamma_2` do not give third points on one
common smaller circle with common fixed endpoints. They give two separate
three-point problems:

    (delta_1 a, delta_1 b, c_1) on Norm=N/Norm(gamma_1),
    (delta_2 a, delta_2 b, c_2) on Norm=N/Norm(gamma_2).

The Markov chain has no common descended space in which these states can mix,
compete, or consume one shared geometric budget.

## 3. Critical-scale invariance

In the balanced squarefree regime one has

    Norm(g)=N^(1/2+o(1)),
    Norm(gamma)=N^(1/4+o(1)),
    Norm(delta)=N^(1/4+o(1)).

Hence

    R_gamma=N^(3/8+o(1)),
    L_gamma=N^(1/8+o(1))=R_gamma^(1/3+o(1)).

So every state lands exactly at the classical Jarnik critical scale. This
explains the bounded multiplicity of a single state, but it supplies no
interaction among distinct states.

## 4. Equivalence with a star bound

For a fixed nearly critical edge `(i,j)`, every other cluster point `k`
produces one state `gamma_ijk`. The critical-scale three-point bound gives only
bounded multiplicity per state. Therefore

    number of continuations from (i,j)
      is comparable, up to O_C(1),
    to number of realizable states gamma_ijk.

Consequently a uniform bound on the number of realizable states for every
nearly critical edge immediately gives the uniform endpoint theorem after the
balanced-clique extraction.

Conversely, the uniform endpoint theorem trivially bounds the number of such
states. Thus, modulo the already-established balanced-edge and bounded
state-multiplicity reductions, the fixed-edge divisor-realizability theorem is
Equivalent to the original uniform-bound problem.

It is not an independent lower-dimensional theorem.

## 5. Orientation-code interpretation

In the squarefree model, let `A_ij` be the set of conductor primes on which
`z_i,z_j` use the same orientation. For a continuation `k`, the triple gcd is
the oriented product over those primes in `A_ij` where `k` uses that common
orientation.

Thus

    k -> gamma_ijk

is exactly the projection of the third orientation word onto the coordinates
where the edge endpoints agree.

Near-injectivity of this projection is possible in abstract balanced codes.
The Gaussian geometry only says that each projected word, after its own
state-dependent descent, solves a critical three-point problem. Since the
descents are different, the projection words are not coupled by one smaller
configuration.

## 6. Concrete obstruction to the Markov program

The proposed stochastic mechanism required one of the following:

1. nested divisor states along overlapping edges;
2. a common descended circle for different states;
3. an amortized loss of entropy when the state changes;
4. a shared geometric or arithmetic budget consumed by distinct states.

None is present:

- triple gcds on overlapping edges are not nested;
- the descended circle and endpoints depend on the state;
- state entropy can be as large as continuation entropy;
- each state pays its own critical-scale three-point cost independently.

Therefore the Markov/divisor-state route does not create an iterative
contraction. It repackages the large endpoint cluster as a large list of
state-dependent critical three-point configurations.

## 7. Verdict

The concrete obstruction is:

> The state normalization is exactly scale invariant at the endpoint, and
> distinct states decouple after descent.

This prevents a Markov contraction or renewal argument. Any theorem uniformly
bounding the realizable states must use a new cross-state identity not present
in the current graph, curvature, gcd, or entropy data. Without such an
identity, the fixed-edge theorem is effectively the original conjecture in
star form.
