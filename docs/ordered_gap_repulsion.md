# Ordered gap repulsion on a lattice circle

## Status

This note develops the ordered-gap viewpoint for lattice points on one short circular arc. It proves an exact local repulsion law for consecutive gaps and shows precisely what it does and does not imply for the endpoint problem.

Let the circle have radius `R`, and let

    z_0, z_1, ..., z_m

be lattice points in positive angular order on an arc of length at most

    C sqrt(R).

Write

    g_i = arg(z_{i+1})-arg(z_i) > 0,
    ell_i = |z_{i+1}-z_i|.

## 1. Exact adjacent-chord determinant

Let

    v_i = z_{i+1}-z_i in Z^2.

The direction of `v_i` is the tangent direction at the midpoint angle of the corresponding chord. Hence the turning angle from `v_i` to `v_{i+1}` is

    (g_i+g_{i+1})/2.

Therefore

    det(v_i,v_{i+1})
      = ell_i ell_{i+1} sin((g_i+g_{i+1})/2).

The determinant is a positive integer, so

    1 <= ell_i ell_{i+1} sin((g_i+g_{i+1})/2).

Using

    ell_i = 2 R sin(g_i/2)

and the short-arc assumption, one obtains for sufficiently large `R`

    ell_i ell_{i+1}(ell_i+ell_{i+1}) >= c_C R,

with an absolute positive constant after fixing the arc constant `C`.

This is the ordered repulsion law.

## 2. Consequences

If `ell_i <= ell_{i+1}`, then

    ell_{i+1} >= c sqrt(R/ell_i).

Thus a bounded physical gap forces each adjacent gap to have order at least `sqrt(R)`. In particular, only `O_C(1)` consecutive gaps can remain bounded as `R` grows.

More generally, writing

    ell_i = R^{alpha_i},

at the exponent level the inequality gives

    alpha_i + alpha_{i+1} + max(alpha_i,alpha_{i+1}) >= 1-o(1).

If `alpha_i <= alpha_{i+1}`, then

    alpha_{i+1} >= (1-alpha_i)/2-o(1).

The fixed point is `alpha=1/3`. Hence gaps below the Jarnik scale `R^{1/3}` force larger neighbours, whereas a chain of gaps all of order `R^{1/3}` is not excluded by this argument.

## 3. Multiscale counting

Let `N_u` be the number of gaps with

    ell_i <= u.

Every such gap must have an adjacent gap at least

    c sqrt(R/u).

After a bounded-overlap charging of small gaps to adjacent large gaps,

    N_u <= C_C sqrt(u)

for `1 <= u <= R^{1/3}`, up to harmless endpoint constants. Indeed the total arc-length budget is `O_C(sqrt(R))`, while every charged large neighbour costs `Omega(sqrt(R/u))`.

This proves a uniform bound for gaps of bounded size and a square-root tail for sub-Jarnik gaps.

At the transition scale `u=R^{1/3}`, it gives only

    N_u = O_C(R^{1/6}),

which matches the classical convex-lattice scale and is not uniform.

## 4. Interaction with the Gaussian residual factorization

For a consecutive pair in the squarefree conductor model, write

    ell_i^2 = (N/P_i) r_i,

where `P_i` is the product of switched rational primes and `r_i` is the norm of the residual Gaussian chord after removing the common conductor factor.

The equal-norm condition implies that `r_i` is even. Indeed if two Gaussian integers have equal norm, the squared norm of their difference is even. Thus the formal residual minimum is `2`, not `1`.

The adjacent determinant law says that two residual-minimal or bounded-residual gaps cannot occur consecutively unless one of the conductor-distance excesses is large enough to make the neighbouring physical gap macroscopic.

This confirms the qualitative intuition that the critical lower spacing cannot be attained repeatedly along the arc.

## 5. Why this does not yet prove the endpoint bound

The repulsion law controls very small gaps, but it has a stable scale at

    ell_i approximately R^{1/3}.

A chain of `O(R^{1/6})` gaps at that scale is compatible with the determinant and total-length inequalities. To obtain a uniform bound, one needs an arithmetic theorem excluding a long run of near-`R^{1/3}` consecutive chords on one fixed circle.

Equivalently, the missing ordered statement is stronger than adjacent repulsion:

> For a long ordered chain of lattice points on one endpoint arc, a positive proportion of consecutive triples must have determinant substantially larger than `1`, in a way whose total gain grows with the chain length.

The determinant identity alone only gives the lower bound `1` independently at each turn.

## 6. Revised target

The ordered-gap approach reduces the conjecture to an integer-curvature accumulation problem.

Define

    K_i = det(z_{i+1}-z_i, z_{i+2}-z_{i+1}) in Z_{>0}.

Then

    K_i = ell_i ell_{i+1} sin((g_i+g_{i+1})/2).

A sufficient theorem would be a superlinear lower bound

    sum_i log K_i >= c m log m - O_C(m),

or any comparable statement forcing cumulative integer curvature to grow faster than linearly in the number of points. Combined with the total arc-length and curvature budget, this would yield a uniform bound.

At present this curvature accumulation theorem is unproved. The contribution of the ordered approach is a genuine local repulsion theorem and a precise identification of the remaining stable scale `R^{1/3}`.
