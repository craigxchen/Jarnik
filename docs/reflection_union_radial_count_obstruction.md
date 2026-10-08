# Full reflection unions have a cardinality-dependent endpoint cost

A whole-set reflection union cannot amplify an unbounded point count while
keeping a bounded endpoint constant. After the least integral realization,
its reflection axis is a coordinate axis or diagonal. Radial lattice spacing
then gives a uniform count for the symmetric union itself. This elementary
obstruction applies to arbitrary nested prime allocations and arbitrary
Gaussian units, without a balanced-profile hypothesis.

The exact denominator and primitive-union facts below are already proved in
[midpoint_reflection_scaling_obstruction.md](midpoint_reflection_scaling_obstruction.md)
and [rational_reflection_route.md](rational_reflection_route.md). The
additional point here is the radial count and its consequent lower bound on
the cost of every pair-reflection union. It does not bound the original
unsymmetrized cluster.

## 1. Primitive symmetric configurations have only unit reflection axes

Let \(U\subset\mathbb Z[i]\) have common nonzero modulus \(R_U\), Gaussian
gcd one, and be invariant under a rational reflection
\(z\mapsto\rho\bar z\), \(|\rho|=1\). Bezout for the conjugate tuple
shows that \(\rho\) is a Gaussian integer. Its modulus is one, so
\(\rho\in\{1,-1,i,-i\}\). Thus its axis is a coordinate axis or diagonal.

The common squared norm \(R_U^2\) is odd: if it were even, every point
would be divisible by \(1+i\), contradicting the common gcd hypothesis.
Projection onto a coordinate axis has integer values. Projection onto a
diagonal is \((x+y)/\sqrt2\) or \((x-y)/\sqrt2\); because the squared
norm is odd, these integer numerators are all odd, so distinct diagonal
projection levels are separated by at least \(\sqrt2\). In every case the
spacing between distinct radial projection levels is at least one.

Suppose \(U\) consists of \(M\) distinct points in its symmetric containing
arc of angular width \(\Theta<\pi\), centered on the relevant axis ray.
Each radial projection level contains at most two circle points. Thus there
are at least \(\lceil M/2\rceil\) distinct levels and
\[
 R_U(1-\cos(\Theta/2))\ge\lceil M/2\rceil-1.
\]
Writing its actual normalized endpoint constant as
\(C_U=\Theta\sqrt{R_U}\), the inequality \(1-\cos x\le x^2/2\) gives
\[
 \boxed{C_U^2\ge8(\lceil M/2\rceil-1).} \tag{1}
\]
Equivalently \(M\le C_U^2/4+2\). The estimate concerns actual integral
symmetric configurations; no arbitrary rotation of the ambient lattice is
being treated as integral.

## 2. Exact least-radius cost of the pair-reflection union

Let \(Z=\{z_1,\ldots,z_m\}\) be a primitive equal-radius cluster of
squared radius \(N=R^2\). Reflection about the midpoint phase of an anchor
pair \(a,b\) has coefficient
\[
 \rho_{ab}=z_a z_b/N=A/D,
 \qquad (A,D)_{\mathbb Z[i]}=1.
\]
The least integral realization of the full union is
\[
 U=DZ\ \cup\ A\bar Z.
\]
It is primitive, and its least squared radius is exactly
\[
 N_U=NQ_{ab},\qquad Q_{ab}=N(D). \tag{2}
\]
Indeed any common multiplier making the original primitive tuple integral
is Gaussian integral by Bezout; integrality of the reflected tuple then
forces divisibility by \(D\). Conversely \(D\) works, and the resulting
union has gcd one since a common divisor would divide both \(A\) and \(D\).
Its reflection coefficient is \(A/\bar D\), a Gaussian unit by Section 1.

This cost has the exact prime-power expression requested in the proposed
symmetrization route. For source half-angle signed valuations \(t_{i,p}\),
put \(u_p=\min_i t_{i,p}\), \(v_p=\max_i t_{i,p}\), and \(W_p=v_p-u_p\).
The reflected row valuations are \(t_{a,p}+t_{b,p}-t_{i,p}\), hence the
union width is
\[
 W_p+|t_{a,p}+t_{b,p}-u_p-v_p|.
\]
Consequently
\[
 \boxed{Q_{ab}=\prod_{p\equiv1\pmod4}
 p^{|t_{a,p}+t_{b,p}-u_p-v_p|}.} \tag{3}
\]
To check this directly in a primitive source realization, use allocations
\(e_{i,p}\in[0,W_p]\) whose extrema are zero and \(W_p\). The reflection
coefficient has opposite Gaussian valuations
\(e_{a,p}+e_{b,p}-W_p\) and its negative. Its reduced denominator norm
therefore has the absolute exponent in (3). Arbitrary integer depths are
retained. Gaussian units contribute no denominator norm, and there is no
ramified or inert denominator contribution.

## 3. Uniform cost over every choice of anchor pair

Assume \(m\ge3\), and the source lies in an arc of actual angular width
\(\Delta\) with
\(C=\Delta\sqrt R\le2\). The source satisfies the sharp pair-product
rigidity theorem from
[multiplicative_rectangle_separation.md](multiplicative_rectangle_separation.md),
whose forbidden threshold is \(2\sqrt2\).
Thus its intersection with its pair-midpoint reflection consists of exactly
the two anchors: another intersection would give a different unordered
pair with product \(z_a z_b\), including the possible repeated-factor case.
Therefore
\[
 |U|=2m-2. \tag{4}
\]

If \(m\ge3\), the source cannot have primitive squared radius one; hence
\(N\ge5\), and \(\Delta\le2/5^{1/4}<\pi/2\). Let \(L\) be the maximum
angular distance of a source point from the chosen midpoint axis. Since the
axis lies between the source endpoints, \(L\le\Delta\). The full union has
the symmetric angular interval \([-L,L]\), of width at most \(2\Delta<\pi\).
It is the shortest containing arc. From (1) and (4),
\[
 \boxed{C_U\ge\sqrt{8(m-2)}.} \tag{5}
\]
On the other hand, (2) gives
\[
 C_U=2L\sqrt{R_U}\le2C Q_{ab}^{1/4}.
\]
Combining these inequalities proves the completely general cost bound
\[
 \boxed{Q_{ab}\ge\frac{4(m-2)^2}{C^4}.} \tag{6}
\]
If the anchors are the two outermost angular points, \(L=\Delta/2\), so
the union has exactly the original angular width. In that case
\[
 C_U=C Q_{ab}^{1/4},\qquad
 \boxed{Q_{ab}\ge\frac{64(m-2)^2}{C^4}.} \tag{7}
\]
Using any prescribed upper endpoint constant in place of the actual \(C\)
weakens these bounds and remains valid.

## 4. Why averaging, rotation, and nested full symmetrizations do not evade it

Inequality (5) is uniform over anchor pairs. Averaging or minimizing the
primewise cost cannot supply a bounded-constant full union when \(m\)
grows. A common rotation changes neither the least radius of the phase
configuration nor its angular width; once the rotated union is placed in
its least primitive integer realization, its reflection is again a Gaussian
unit and (1) applies.

After any final operation of the form \(V\mapsto V\cup\sigma(V)\), the
result is invariant under that last reflection. Therefore (1) applies
regardless of earlier nested symmetrizations or their prime-power history.
Keeping an increasing full point set forces its final endpoint constant to
grow at least as the square root of its cardinality.

Discarding points after symmetrization can destroy this symmetry and avoids
the radial bound, but then the retained number, its new primitive radius,
and its arc width must all be charged. No count-amplification statement for
that different operation is established here. This is an obstruction to
full-union symmetrization as a proposed proof operation, not a contradiction
to the original short-arc cluster or a solution of the uniform count goal.

Independent exact checks used all six pair reflections of each actual Pell
cluster at indices \(1,11,21\). All 18 unions had exactly six points, common
Gaussian gcd one, a unit reflection coefficient, and exactly three radial
projection levels. Diagonal levels had the required odd numerator parity.
The unrestricted count and cost bounds follow from the proofs above.

The persistent exact Pell verification is
[check_reflection_union_radial_count.py](check_reflection_union_radial_count.py).
