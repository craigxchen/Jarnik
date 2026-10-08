# Fixed-dimensional clique audit for t = 3, 4, 5

## Status

This note tests the fixed-dimensional near-square clique reduction at the first
three nontrivial sizes.  The conclusion is concrete:

* the rational tangent cocycle is exactly parametrized by vertex slopes and adds
  no new rigidity at any fixed t;
* the critical gcd-pattern weights are unique for t=3, form a one-parameter
  parity family for t=4, and form a five-parameter family for t=5;
* therefore increasing t does not make the exponent/gcd system finite or
  zero-dimensional;
* any successful fixed-t argument must use the actual small integer residuals
  in the Gaussian equations, not only norms, gcds, cocycles, or logarithmic
  pattern weights.

The uniform endpoint theorem is not proved here.

---

## 1. Cocycle elimination is tautological

Let the selected points be indexed by 0,...,t-1.  For every oriented edge
choose a nonzero Gaussian integer

    B_ij = X_ij + i Y_ij

such that

    B_ji = conjugate(B_ij)

and suppose

    B_ij B_jk conjugate(B_ik) is real and nonzero

for every triple.  Equivalently, the norm-one phases

    u_ij = B_ij / conjugate(B_ij)

satisfy

    u_ij u_jk = u_ik.

Choose u_i = u_i0.  Then necessarily

    u_ij = u_i / u_j.

Thus every complete cocycle comes from t vertex phases.  Writing a rational
half-angle representative

    A_i = x_i + i y_i

for u_i gives, up to a nonzero real scalar on each edge,

    B_ij ~ A_i conjugate(A_j).

Consequently

    X_ij ~ x_i x_j + y_i y_j,
    Y_ij ~ y_i x_j - x_i y_j.

The Y_ij are the 2 by 2 minors of the 2 by t matrix with columns (x_i,y_i).
All polynomial identities among them are the rank-two Plucker identities.
There is no new algebraic collapse at t=4 or t=5: the cocycle variety remains
positive-dimensional for every t.

In particular:

* t=3: one tangent-addition identity;
* t=4: the single quadratic Plucker relation on each quadruple;
* t=5: the five Plucker relations cutting out Gr(2,5), still of dimension 7
  projectively before the additional norm conditions.

Therefore eliminating the X_ij from the cubic cocycle equations cannot by
itself produce a finite classification.

---

## 2. Prime-pattern coordinates

Root at point 0.  Every squarefree conductor prime has a membership pattern

    s in F_2^(t-1),

recording which of points 1,...,t-1 reverse its Gaussian orientation relative
to point 0.  Let n_s be the product of the rational primes with pattern s and
put

    w_s = log(n_s) / log(N).

Then

    w_s >= 0,
    sum_s w_s = 1.

The normalized logarithmic distance between points i and j is the probability,
under the distribution w, that the corresponding bits differ.  The balanced
clique condition says

    P(s_i != s_j) = 1/2 + o(1)

for every pair, with the root bit fixed to zero.

At exact criticality this says that the random signs

    epsilon_i = (-1)^(s_i)

are unbiased and pairwise orthogonal:

    E epsilon_i = 0,
    E epsilon_i epsilon_j = 0  (i != j).

Hence the critical pattern distributions are exactly the pairwise-independent
distributions on F_2^(t-1).  Fourier analysis on the cube gives the normal
form

    w(s) = 2^{-(t-1)} sum_{T subset {1,...,t-1}} c_T chi_T(s),

where

    c_empty = 1,
    c_T = 0 for |T|=1 or 2,

and only Fourier coefficients of degree at least 3 remain free, subject to
positivity.

This gives exact dimensions for small t.

---

## 3. t = 3

There are two non-root bits and four pattern blocks

    n_00, n_10, n_01, n_11.

The three critical distances are

    log d_1       = log n_10 + log n_11 = H/2,
    log d_2       = log n_01 + log n_11 = H/2,
    log q_12      = log n_10 + log n_01 = H/2.

Together with total mass H, these equations have the unique solution

    log n_00 = log n_10 = log n_01 = log n_11 = H/4.

Thus a critical triple forces four conductor blocks of quarter-height.
This is exactly the balanced three-block factorization already obtained by the
Gaussian residual method.

After removal of forced common factors, one gets three near-real Gaussian
integers with norms

    n_10 n_01,
    n_10 n_11,
    n_01 n_11,

all of size N^(1/2+o(1)), and one primitive three-term Gaussian relation with
subpower residual coefficients.

This is a genuine simplification, but it is not finite: the four quarter-size
squarefree blocks can vary without bound.  Primitive three-term Gaussian
relations with variable pairwise-coprime supports are flexible.

Verdict for t=3:

    exact normal form, but not enough rigidity.

---

## 4. t = 4

There are three non-root bits and eight pattern blocks.  Exact criticality
forces all Fourier coefficients of degrees one and two to vanish.  The only
remaining coefficient is the parity coefficient tau = c_{123}.  Hence

    w(s) = (1 + tau (-1)^(s_1+s_2+s_3)) / 8,
    -1 <= tau <= 1.

Equivalently, every even-parity pattern has weight

    (1+tau)/8,

and every odd-parity pattern has weight

    (1-tau)/8.

So t=4 does not produce a finite list.  It produces a continuous one-parameter
family.

The endpoints tau=+1 and tau=-1 are the two parity-code distributions supported
on four patterns.  The midpoint tau=0 is the uniform distribution on all eight
patterns.  Every value in between remains exactly pairwise critical.

The four-point cocycle/Plucker equation is automatically satisfied by four
vertex half-angle slopes.  At the valuation level it only says that the three
opposite-edge products in Ptolemy are tropically balanced; the parity family
already has that property.

Verdict for t=4:

    one extra point exposes parity, but leaves a full one-parameter critical
    family.  No zero-dimensional classification appears.

---

## 5. t = 5

There are four non-root bits and sixteen pattern blocks.  Exact criticality
kills Fourier coefficients of degrees one and two.  The free coefficients are

    c_123, c_124, c_134, c_234, c_1234.

Thus the affine space of exact critical signed distributions has dimension
five before imposing nonnegativity.  The feasible set contains a neighborhood
of the uniform distribution, so it is genuinely five-dimensional.

This is already worse than t=4.  The ten pairwise distance equations do not
come close to determining the sixteen pattern masses.

The ten Gaussian edge variables are also not independent: they are generated
by five vertex phases.  Their imaginary parts are Plucker coordinates of a
rank-two 2 by 5 matrix.  The resulting Grassmannian relations organize the
system but do not reduce it to finitely many components.

Verdict for t=5:

    the exponent/gcd geometry becomes less rigid; five higher-order correlation
    parameters survive.

---

## 6. Why larger fixed t does not repair this

For t points, exact pairwise criticality leaves every Fourier coefficient of
pattern degree at least three unconstrained.  Their number is

    2^(t-1) - 1 - (t-1) - binom(t-1,2),

which grows exponentially with t.

Thus passing to larger fixed cliques cannot simplify the prime-pattern weights
using pairwise near-criticality alone.  It creates more higher-order pattern
freedom than equations.

Similarly, the cocycle equations always reduce to t vertex phases.  They do not
become overdetermined as the number of edges grows, because all edge relations
are syzygies of the same rank-two parametrization.

---

## 7. What information is genuinely left

After this audit, the fixed-dimensional problem has only one possible source of
additional rigidity: the actual small integer residuals.

For each edge,

    B_ij = X_ij + i Y_ij,
    |Y_ij| = N^o(1),
    Norm(B_ij) = q_ij = N^(1/2+o(1)).

The logarithmic pattern equations see only Norm(B_ij).  The cocycle sees only
the edge phases.  What neither abstraction uses is that the same small integers
Y_ij must simultaneously satisfy

    Y_ij ~ det((x_i,y_i),(x_j,y_j))

with exact divisibility by the pattern blocks and exact Gaussian factorization
on every overlapping triple and quadruple.

A tractable fixed-t theorem would therefore have to bound primitive integral
rank-two configurations with:

1. all Plucker coordinates subpower;
2. all column norms prescribed squarefree divisors of one N;
3. every pairwise symmetric-difference norm near sqrt(N);
4. the forced local prime factors removed from each Plucker coordinate.

That is an integral small-Plucker-coordinate problem, not a cocycle-elimination
problem.

---

## 8. Final verdict

Trying t=3,4,5 gives a clear answer:

    t=3  -> unique quarter-block profile, flexible Gaussian S-unit equation;
    t=4  -> one-parameter parity profile;
    t=5  -> five-parameter pairwise-independent profile.

Therefore the proposed strategy "increase t and eliminate the cocycle" does not
simplify the problem into a finite algebraic classification.

The viable fixed-dimensional reformulation is narrower:

> Prove that a primitive integral 2 by t matrix cannot have all forced-reduced
> Plucker coordinates of size N^o(1) while its column norms and pairwise
> symmetric-difference norms all occupy the critical N^(1/2+o(1)) scale.

The smallest value worth attacking is t=4, because t=3 is a flexible
three-term equation and t=4 is the first case with a nontrivial Plucker relation.
But the parity parameter tau shows that one must use the exact residual
integers, not just their sizes or valuations.
