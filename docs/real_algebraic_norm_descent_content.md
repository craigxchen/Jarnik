# Primitive content of a real-algebraic norm descent

This note concerns the certified four-row real-algebraic polynomial profile
of [four_row_real_polynomial_profile_certificate.md](four_row_real_polynomial_profile_certificate.md).
It gives a uniform Gaussian-content calculation for **rational parameters**.
The calculation does not establish that the certified roots have distinct
Galois orbits, nor does it turn that profile into integer points on short
arcs.

## The five actual circle polynomials

Let \(K\subset\mathbb R\) contain the real and imaginary coordinates of
the fifteen certified roots \(z_T\), let \(L=K(i)\), and put
\(d=[L:\mathbb Q(i)]\). For nonempty \(T\subseteq[4]\), set

\[
 J_T(t)=N_{L/\mathbb Q(i)}(t-z_T)\in\mathbb Q(i)[t].
\]

Each \(J_T\) is monic of degree \(d\). Choose fixed positive integers \(M_T\)
so that the homogeneous forms

\[
 B_T(A,B)=M_T B^d J_T(A/B)\in\mathbb Z[i][A,B]
\]

have integral coefficients. Write a bar for Gaussian coefficient
conjugation, and define

\[
 Z_0=\prod_T B_T,\qquad
 Z_i=\prod_{T\not\ni i}B_T\prod_{T\ni i}\overline{B_T}
 \quad(1\leq i\leq4).                                      \tag{1}
\]

At every real parameter pair \((a,b)\), the five values have the same
modulus. The four row products \(Q_i=\prod_{T\ni i}B_T\) are **not** the
five circle points: in particular \(B_{[4]}\) divides every \(Q_i\).
Their gcd therefore cannot measure primitive content of (1).

## Exact content, including correction primes

Factor the thirty oriented forms \(B_T,\overline{B_T}\) over the Gaussian
UFD into primitive irreducible homogeneous forms \(\Phi_O\), indexed by
the \(\mathbb Q(i)\)-Galois orbits \(O\) of their roots. The field norm of
\(t-z_T\) is a power of the minimal polynomial of \(z_T\), so each
oriented label has just one such orbit, possibly with multiplicity. Write

\[
 Z_j=c_j\prod_O\Phi_O^{e_{jO}},\quad
 \mu_O=\min_{0\leq j\leq4}e_{jO},\quad
 D=\prod_O\Phi_O^{\mu_O},                                  \tag{2}
\]

where \(c_j\in\mathbb Z[i]\setminus\{0\}\). Distinct \(\Phi_O\) have
nonzero homogeneous resultants \(R_{OP}\in\mathbb Z[i]\). Let
\(E=\max_{j,O}(e_{jO}-\mu_O)\), and take

\[
 C=\Bigl(\prod_{j=0}^4c_j\Bigr)
   \prod_{O<P}R_{OP}^{E}\ne0.                               \tag{3}
\]

For every primitive pair \((a,b)\in\mathbb Z^2\) for which the values in
(1) are nonzero, their Gaussian gcd \(G(a,b)\) satisfies, up to a unit,

\[
       D(a,b)\mid G(a,b)\mid C D(a,b).                         \tag{4}
\]

This is a prime-power statement, not just a bound on the set of prime
divisors. To prove it, fix a Gaussian prime \(\pi\), put
\(x_O=v_\pi(\Phi_O(a,b))\), and choose \(O_0\) with maximal \(x_O\).
The homogeneous resultant identity gives

\[
 \min(x_O,x_P)\leq v_\pi(R_{OP})\quad(O\ne P),                \tag{5}
\]

because at least one of \(a,b\) is a \(\pi\)-adic unit. Choose a row \(j\)
for which \(e_{jO_0}=\mu_{O_0}\). After subtracting the valuation of \(D\),
the valuation of \(Z_j\) is at most

\[
 v_\pi(c_j)+E\sum_{P\ne O_0}x_P
 \leq v_\pi(C).
\]

The first divisibility in (4) follows directly from (2). This proof
charges primes in leading coefficients, denominators, and collisions at
specialized parameters to the same fixed \(C\); it does not discard them.
Consequently \(G/D\) has only finitely many possible values up to a unit.

If the thirty oriented roots have distinct \(\mathbb Q(i)\)-Galois
orbits, every \(\mu_O=0\): the anchor \(Z_0\) omits each conjugate form,
and an incident row \(Z_i\) omits each unconjugated form. Hence

\[
                         |G(a,b)|\leq |C|.                   \tag{6}
\]

Separation of the thirty roots at the certified real embedding does not
prove this orbit hypothesis. Two separated roots can be Galois conjugate.

## What a collision can remove

Coefficient conjugation pairs orbit classes \(O,\bar O\). If
\(O=\bar O\), the orbit contains both orientations of every label it
contains. Each point in (1) then takes one of those two equal factors
per label, so \(e_{0O}=\cdots=e_{4O}\): this factor is entirely in \(D\)
and supplies no primitive core.

For \(O\ne\bar O\), let \(m_O\) be the number of oriented labels in \(O\),
with field-norm multiplicities included. Every point takes exactly one
orientation from each underlying label, so

\[
 e_{jO}+e_{j\bar O}=m_O\qquad(0\leq j\leq4).          \tag{7}
\]

At each Gaussian prime \(\pi\), select an orbit \(O_\pi\) maximizing
\(x_O=v_\pi(\Phi_O(a,b))\). Equation (5) bounds every other orbit depth
by a fixed resultant valuation, even when \(\pi\) divides that
resultant. The anchored row-valuation vector, unchanged by dividing
all five points by their common gcd, is

\[
 v_{O_\pi}k_\pi+\epsilon_\pi,\quad
 v_O=(e_{1O}-e_{0O},\ldots,e_{4O}-e_{0O}),\quad
 k_\pi=x_{O_\pi},                                    \tag{8}
\]

where the error \(\epsilon_\pi\) contains nondominant orbit depths
and the fixed coefficient valuations. More precisely, for a fixed
constant depending on (2),

\[
 \sum_\pi\|\epsilon_\pi\|_\infty\log|\pi|
 \ll 1+\sum_j\log|c_j|+\sum_{O<P}\log|R_{OP}|=O(1).
                                                        \tag{9}
\]

The same estimate shows that assigning each Gaussian prime to its
maximal orbit loses only \(O(1)\) total logarithmic depth from every
orbit. This retains arbitrarily high powers at fixed resultant primes.
An inert prime cannot carry an unbounded non-self-conjugate depth:
it would divide \(\Phi_O\) and \(\Phi_{\bar O}\) to the same depth,
contradicting their fixed resultant. Conjugate split primes have
opposite anchored vectors by (7).

Define the correction budget of a proposed Boolean profile to be the
sum, over Gaussian primes, of the maximum-norm difference between the
actual anchored valuation vector and the profile's vector, weighted
by \(\log|\pi|\). A single coherent core on a nonempty cut \(T\) has
vector \(s k'\mathbf1_T\) at a prime, where \(s\in\{-1,1\}\) and
\(k'\geq0\). Distinct cores have disjoint rational-prime supports.
Multiplying the five points by a fixed number of Gaussian rational
correction factors whose combined numerator and denominator logarithmic
heights are \(o(\log a)\) gives precisely an \(o(\log a)\) budget in
this sense.
If \(v_O\ne0\) and is not on a Boolean ray, its distance \(\delta_O\)
from the union of those rays is positive. Equation (8) charges at
least \(\delta_Ok_\pi-\|\epsilon_\pi\|_\infty\) depth at every prime
assigned to \(O\). Summing gives a fixed positive fraction of
\(\log N(\Phi_O(a,b))\), up to \(O(1)\). If \(v_O=0\), the orbit pair
supplies no core. If \(v_O\) is on a Boolean ray, prime mass assigned
to any other cut is \(o(\log a)\) under the stated correction budget,
even when the nonzero entries have magnitude greater than one.
Splitting their prime-power depths cannot create disjoint supports.

Therefore, along \(a\to+\infty\) with fixed \(b\ne0\), a full set of
fifteen **distinct, prime-disjoint, balanced Boolean cores** cannot
survive an orbit collision with \(o(\log a)\) correction budget.
Each non-self-conjugate orbit pair has positive logarithmic mass
\(2\deg(\Phi_O)\log a+O(1)\), supplies at most one cut under that
budget. There are thirty oriented labels; a non-self-conjugate pair
uses at least two of them, and a self-conjugate class supplies no cut.
Any collision therefore leaves fewer than fifteen active pairs.
This conclusion uses the stated full-core and correction hypotheses.
It does not rule out a different, smaller cut profile after large
systematic cancellation.

## Exact chord cost after primitive reduction

For \(b=1\) and \(t=a\to+\infty\), put
\(A_i=\prod_{T\not\ni i}B_T(t,1)\) and
\(Q_i=\prod_{T\ni i}B_T(t,1)\). Suppose
\(q_i=\deg\operatorname{Im}Q_i\) is finite. Equation (1) gives

\[
 Z_i-Z_0=A_i(\bar Q_i-Q_i),\qquad
 |Z_i-Z_0|=2|A_i|\,|\operatorname{Im}Q_i|.            \tag{10}
\]

Write \(g=\deg D\). Since \(|A_i|\asymp t^{7d}\),
\(|Z_0|\asymp t^{15d}\), and (4) gives
\(|G|\asymp t^g\), the primitive chord and radius
\(R=|Z_0/G|\) obey

\[
 \frac{|(Z_i-Z_0)/G|}{\sqrt R}
       \asymp t^{q_i-(d+g)/2}.                        \tag{11}
\]

Thus bounded content cannot rescue a row with \(q_i>d/2\) in the
distinct-orbit case. In that case \(\operatorname{Im}Q_i\) cannot vanish
identically: \(Q_i\) and \(\bar Q_i\) have disjoint irreducible supports,
whereas a real polynomial would equal its conjugate. For the original
row \(P_i=X_i+iY_i\), the
coefficient of \(t^{8d-8}\) in the imaginary part of
\(N_{L/\mathbb Q(i)}P_i\) is
\(\operatorname{Tr}_{K/\mathbb Q}(Y_i)\). If this trace is nonzero,
\(q_i=8d-8\); for \(d\geq2\), (11) diverges. A zero trace requires
checking lower coefficients. The borderline \(d=2,q_i=1\) is possible
for an individual row: with \(K=\mathbb Q(\sqrt2)\) and
\(P=(t^2+1)t^6+\sqrt2(t-i)\), the imaginary part of
\(N_{K(i)/\mathbb Q(i)}P\) is \(4t\). This single-row example is not a
full fifteen-block profile.

Formula (11) assumes a fixed denominator and a nonzero imaginary
polynomial. It does not exclude cancellation on moving
rational-parameter sequences. It also does not bound the arc by itself:
it gives a necessary chord condition for the five points in (1).

## Algebraic-unit paths

For an algebraic unit \(\theta\in K\), the Gaussian sequence
\(N_{L/\mathbb Q(i)}(\theta^n-z_T)\) is generally **not** the value
\(J_T(\theta^n)\): its conjugate factors contain
\(\sigma(\theta)^n\), not a common indeterminate \(\theta^n\).
The one-variable resultant proof above therefore does not apply
automatically. A fixed progression and a rank-one Laurent
parametrization of all conjugate units would permit a new resultant
calculation, provided the resulting oriented Laurent polynomials have
no common root after their genuine collisions are grouped. For a
higher-dimensional torus of conjugate powers, coprimality of two
functions alone gives no constant in their ideal; simultaneous
high-depth congruences at an intersection require separate analysis.
No bounded primitive-content claim, periodicity exclusion, or
sparse-subsequence exclusion is made for a general algebraic-unit path.
