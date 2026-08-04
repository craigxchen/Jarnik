# Fixed-dimensional near-square clique reduction

## Status

This note replaces the moving prime-pattern formulation by one fixed finite
Diophantine system. It is a reduction, not a proof of the uniform endpoint
bound.

## 1. Rooted squarefree model

Assume the squarefree diffuse case and choose one point z_0 in a balanced
clique of fixed size t. For every other point z_i, let S_i be the set of
rational conductor primes at which z_i uses the opposite Gaussian orientation
from z_0. Let

    d_i = product_{p in S_i} p.

Choose the oriented Gaussian product

    A_i = product_{p in S_i} pi_p = x_i + i y_i,

so Norm(A_i)=d_i. Up to a common Gaussian divisor,

    z_i/z_0 = A_i/conjugate(A_i).

The chord identity gives

    |z_i-z_0|^2 = (N/d_i) |A_i-conjugate(A_i)|^2
                = 4 (N/d_i) y_i^2.

For a balanced clique,

    d_i = N^{1/2+o(1)}.

The endpoint chord bound |z_i-z_0|^2 <= C^2 N^{1/2} therefore implies

    |y_i| = N^{o(1)},
    |x_i| = N^{1/4+o(1)}.

Thus every selected point is represented by a primitive Gaussian divisor of N
whose argument is controlled by a subpower numerator over a quarter-power
denominator.

## 2. Pairwise symmetric-difference divisors

For i,j define

    g_ij = gcd(d_i,d_j),
    q_ij = d_i d_j / g_ij^2.

Since N is squarefree, q_ij is the product of the conductor primes on which
points i and j differ. Balanced clique extraction gives

    q_ij = N^{1/2+o(1)}

for every selected pair.

Remove the common oriented Gaussian factors from A_i and A_j. The remaining
Gaussian product B_ij has

    Norm(B_ij)=q_ij

and the pairwise chord identity gives

    B_ij = X_ij + i Y_ij,
    q_ij = X_ij^2 + Y_ij^2,
    |Y_ij| = N^{o(1)}.

Therefore every symmetric-difference divisor q_ij is subpower-close to a
perfect square:

    q_ij = X_ij^2 + N^{o(1)}.

This holds simultaneously for all pairs in one fixed clique.

## 3. Exact cocycle compatibility

The norm-one quotients satisfy

    (z_i/z_j)(z_j/z_k)=z_i/z_k.

Writing

    z_i/z_j = B_ij/conjugate(B_ij)

shows that

    B_ij B_jk conjugate(B_ik)

is a nonzero real Gaussian integer. Hence the near-square representations are
not independent. If B_ab=X_ab+iY_ab, their real and imaginary parts satisfy an
exact cubic identity, equivalently the tangent-addition law

    Y_ik/X_ik
      = (Y_ij/X_ij + Y_jk/X_jk)
        / (1 - (Y_ij Y_jk)/(X_ij X_jk)),

whenever the displayed denominators are nonzero, after choosing consistent
orientations.

Thus the reduced object is a finite rational tangent cocycle, not an arbitrary
collection of almost squares.

## 4. Finite theorem sufficient for the endpoint bound

It is enough to prove the following statement for one fixed t and one fixed
small epsilon>0.

> There is no arbitrarily large squarefree N admitting squarefree divisors
> d_0=1,d_1,...,d_{t-1} and coherently oriented Gaussian integers B_ij such that
>
> 1. N^{1/2-epsilon} <= q_ij <= N^{1/2+epsilon} for every i!=j;
> 2. q_ij=Norm(B_ij)=X_ij^2+Y_ij^2;
> 3. |Y_ij| <= N^epsilon;
> 4. B_ij B_jk conjugate(B_ik) is real for every triple;
> 5. q_ij=d_i d_j/gcd(d_i,d_j)^2.

Balanced-clique extraction supplies this system with epsilon tending to zero
for every fixed t. A contradiction for one fixed pair (t,epsilon) therefore
proves the uniform endpoint bound.

## 5. Why this is a genuine simplification

The original problem has a moving number of Gaussian prime coordinates and a
moving pattern matrix. The reduced problem has only

    t(t-1)/2

Gaussian variables B_ij and t squarefree divisors d_i, where t is fixed once
and for all. All prime information appears only through gcds among the d_i.

There are no determinant inequalities, exceptional algebraic sets, entropy
increments, or growing-rank pattern kernels left.

The remaining problem is a fixed-dimensional Diophantine classification of a
coherent clique of near-real Gaussian integers of norm about sqrt(N).

## 6. Immediate lines of attack

The most direct attacks are now finite:

1. eliminate X_ij using the cubic cocycle identities and derive polynomial
   equations in the subpower variables Y_ij and the divisors q_ij;
2. use q_ij=d_i d_j/gcd(d_i,d_j)^2 to classify the possible gcd matrix of a
   fixed t-tuple;
3. combine the gap between consecutive squares near N^{1/2} with the exact
   tangent cocycle;
4. seek a fixed t for which the resulting algebraic variety has only degenerate
   components, corresponding to repeated points or a heavy common divisor.

This is the active reduced problem.
