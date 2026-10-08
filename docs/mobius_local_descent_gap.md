# A discrete local gap around conformal endpoint maps

## Status

Constrained minimization does not admit an infinitesimal departure from the
conformal maps.  For any primitive rational circle configuration of least
squared radius \(N\) with at least three distinct points, every
nonconformal rational projective image of least squared radius \(N'\le N\)
has normalized operator distance at least
\[
 \boxed{\frac1{8N+1}}
\]
from the conformal and anticonformal locus.  The conclusion remains true
after imposing any target endpoint bound, including \(C'\le2\).

More quantitatively, a nonconformal perturbation of normalized size
\(\varepsilon<1\) satisfies
\[
 \boxed{\displaystyle
 N'\ge\frac{(1-\varepsilon)^2}{64\varepsilon^2N}.}        \tag{1}
\]
Thus a rational shear by \(p/q\) has a quadratic denominator jump:
\[
 N'\ge\frac{(q-|p|)^2}{64p^2N}\qquad(0<|p|<q).           \tag{2}
\]
This rules out an infinitesimal departure from the conformal locus while
keeping the radius nonincreasing. It neither constructs
nor excludes admissible maps outside this shrinking neighborhood.

## 1. Discreteness of bounded-radius rational phases

Normalize one labeled phase of a primitive rational circle configuration to
one.  Every other relative phase can be written
\[
 q=\frac{H}{\bar H},\qquad H=(x,y)\in\mathbb Z[i],
 \qquad \gcd(x,y)=1.
\]
The Gaussian gcd of \(H\) and \(\bar H\) has norm one or two.  The norm of
the reduced denominator of \(q\) is consequently \(N(H)\) or \(N(H)/2\).
Since this denominator norm divides the least squared radius, a
configuration of least squared radius at most \(X\) admits representatives
satisfying
\[
 N(H)\le2X.                                               \tag{3}
\]
This statement includes both Gaussian parity classes.

For two distinct rational phases \(q=H/\bar H\) and
\(q'=K/\bar K\), direct subtraction gives
\[
 |q-q'|
 =\frac{2|\det(H,K)|}{\sqrt{N(H)N(K)}}.                  \tag{4}
\]
The determinant is a nonzero integer.  If the phases belong respectively
to configurations of least radii \(N\) and \(N'\), equations (3)--(4)
therefore give the uniform separation
\[
 \boxed{\displaystyle
 |q-q'|\ge\frac1{\sqrt{NN'}}.}                            \tag{5}
\]
The phases in (5) may come from different configurations.  Only their
individual denominator bounds are used.

## 2. Projective displacement near a conformal map

For a nonsingular real matrix \(M\), define its scale-free distance from the
conformal and anticonformal locus by
\[
 \varepsilon(M)=
 \inf_{\lambda>0,\ R\in O(2)}
 \|\lambda RM-I\|_{\rm op}.                              \tag{6}
\]
Left multiplication by \(\lambda R\) changes the output only by a common
rotation or reflection and a scalar.  It preserves the least target radius,
all relative phase distances, and the endpoint constant.

This distance has an exact condition-number formula.  If
\(\sigma_1\ge\sigma_2>0\) are the singular values and
\(\kappa=\sigma_1/\sigma_2\), then
\[
 \boxed{\displaystyle
 \varepsilon(M)
 =\frac{\sigma_1-\sigma_2}{\sigma_1+\sigma_2}
 =\frac{\kappa-1}{\kappa+1}.}                            \tag{6a}
\]
Equivalently, for \(T=\operatorname {tr}(M^TM)\) and
\(\delta=\det M\),
\[
 \varepsilon(M)^2=\frac{T-2|\delta|}{T+2|\delta|}.        \tag{6b}
\]
Indeed, choose the left orthogonal factor so that the resulting matrix is
orthogonally diagonalizable with eigenvalues \(\sigma_1,\sigma_2\).  The
optimal scale is \(\lambda=2/(\sigma_1+\sigma_2)\), which makes the two
spectral errors equal.  Conversely, singular-value perturbation gives
\[
 \|\lambda RM-I\|_{\rm op}
 \ge\max_i|\lambda\sigma_i-1|
\]
for every \(\lambda,R\), proving optimality.

Fix a representative
\[
 A=I+E,\qquad \|E\|_{\rm op}\le\varepsilon<1.
\]
For every unit real direction \(u\), put \(v=Au/\|Au\|\).  Since
\(\|Au\|\ge1-\varepsilon\),
\[
 \|v-u\|
 \le\frac{2\varepsilon}{1-\varepsilon}.                  \tag{7}
\]
Passing from a direction to its unit-circle phase doubles its angle and
gives
\[
 |q(Au)-q(u)|
 \le\frac{4\varepsilon}{1-\varepsilon}.                  \tag{8}
\]

Now normalize both source and target phases by their labeled anchor.  The
anchor itself incurs the same error as every other row, so the triangle
inequality in (8) gives
\[
 \left|
 \frac{q(Au_i)}{q(Au_0)}
 -\frac{q(u_i)}{q(u_0)}
 \right|
 \le\frac{8\varepsilon}{1-\varepsilon}.                  \tag{9}
\]
Although \(R\) in (6) need not be rational, it cancels from these relative
phases, up to conjugation in the reflection case.  Hence the target phases
in (9) retain their rational representatives and their least-radius bound.

## 3. The radius jump

Suppose \(M\) is nonconformal and the source contains at least three
distinct projective directions.  The normalized relative phases in (9)
cannot all equal their source values.  Otherwise a common output rotation
or reflection would make the projective action fix three distinct real
lines, forcing that action to be the identity.  The original map would then
be conformal or anticonformal.

Choose a label where the two relative phases differ.  Combining the lower
bound (5) with (9) gives
\[
 \frac1{\sqrt{NN'}}
 \le\frac{8\varepsilon}{1-\varepsilon}.
\]
Rearrangement proves (1).  In particular, \(N'\le N\) implies
\[
 1-\varepsilon\le8N\varepsilon,
\]
and hence
\[
 \boxed{\displaystyle
 \varepsilon(M)\ge\frac1{8N+1}.}                         \tag{10}
\]
If the infimum in (6) is not attained, apply the argument with every
\(\varepsilon>\varepsilon(M)\) and pass to the limit.

The target endpoint constraint was not used.  Therefore the conformal locus
is an isolated component of the set
\[
 \{M:N'(M)\le N,\ C'(M)\le2\}
\]
in the normalized operator topology.

There is also a global discreteness statement for the constrained search.
For fixed source labels and fixed \(N\), only finitely many normalized
target relative-phase tuples have least squared radius at most \(N\):
equation (3) leaves only finitely many primitive integer rows.  Three
distinct source directions determine a real projective map from any three
distinct labeled target directions.  Hence, modulo output rotations,
reflections, and scalars,
\[
 \{M:N'(M)\le N\}
\]
is finite.  Adding \(C'(M)\le2\) only removes elements.  Thus a constrained
minimum exists as a finite arithmetic choice, but it need not contain a
nonconformal element.

## 4. Explicit rational shears

Take
\[
 S_{p/q}=
 \begin{pmatrix}1&p/q\\0&1\end{pmatrix},
 \qquad 0<|p|<q,\qquad (p,q)=1.
\]
The displayed representative has
\[
 \|S_{p/q}-I\|_{\rm op}=|p|/q.
\]
It is nonconformal, so (1) immediately becomes (2).  In particular,
\[
 q>(8N+1)|p|\quad\Longrightarrow\quad N'>N.              \tag{11}
\]
Thus making the rational perturbation smaller by increasing its denominator
eventually forces the primitive target radius upward.  Ordinary row
contents do not evade the conclusion: they have already been removed in
the least-radius reconstruction used in (3)--(5).

The \(1/N\) exclusion scale is conservative.  It is sufficient because the
determinant in (4) is integral; no fairness, Gaussian-unit extraction, or
source endpoint lower bound is assumed.

## 5. Relation to four and eight rows

For four points, high-core projective self-maps of Pell clusters give
nonconformal examples with exactly unchanged radius.  Equation (10) says
that these maps must lie outside the local conformal ball; it does not
contradict them.  The product-minimizing critical shear in
[the reduction obstruction](mobius_product_reduction_obstruction.md) also
lies far outside this ball.

For \(m\ge8\), the same local gap holds without losing any source or target
points.  It does not prove that a nonconformal constrained minimizer exists,
because the feasible set may contain only the conformal locus, or may have
additional components at distance at least \(1/(8N+1)\). A descent argument
based on this constrained search still requires global arithmetic
information about those distant components.

When both endpoint constants are at most two, the existing full-unit
condition-number theorem is much stronger.  It gives
\[
 \kappa\ge N^{z_m}/512,\qquad z_m\ge1/66
\]
for every nonconformal radius-nonincreasing map with \(m\ge8\).  By (6a),
whenever \(N^{z_m}>512\),
\[
 \varepsilon(M)\ge
 \frac{N^{z_m}/512-1}{N^{z_m}/512+1}.                    \tag{12}
\]
This tends to one, expressing the necessary approach to rank one.  The
endpoint-free local gap (10) therefore does not tighten the known admissible
critical regime; its role is to rule out a small-variation construction
without any endpoint or point-count hypothesis.

The finite arithmetic checks are in
[the persistent checker](check_mobius_local_descent_gap.py).
