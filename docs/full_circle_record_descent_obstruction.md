# Large full circles can have no nonconformal radius descent

This is an endpoint-free obstruction to constructing a rational projective
descent from cardinality alone. There are primitive integral circle
configurations with arbitrarily many points for which every rational
projective image of nonincreasing least radius is conformal or anticonformal.
They occupy full circles, so they do not satisfy the fixed short-arc
hypothesis as their radii grow. The result does not settle the existence
question for endpoint clusters.

## 1. Conventions and the full-circle symmetry group

Let
\[
 S_N=\{z\in\mathbb Z[i]:\operatorname{Norm}(z)=N\},\qquad
 r_2(N)=|S_N|,
\]
where N is a positive integer represented as a sum of two squares. Choose
\(z_0\in S_N\) and consider the rational circle phases
\[
 Q_N=\{z/z_0:z\in S_N\}.
\]
Represent these by half-angle directions \(q=H/\bar H\), so ordinary
real projective transformations of H act as circle Möbius transformations
of q. Opposite lattice points give distinct phases q and -q; they are not
identified as ordinary radial lines.

The stabilizer \(G_N\) of \(Q_N\) in \(\mathrm{PGL}_2(\mathbb Q)\)
is finite: its action on the finite phase set is faithful since three
distinct half-angle directions determine a projective map. The subgroup
\(G_N^+\) with positive determinant preserves cyclic order. Its action
therefore embeds it in the cyclic group of cyclic shifts of the ordered
phase set, so \(G_N^+\) is cyclic.

The quarter-turn q to iq preserves \(Q_N\) and is represented on H by
\[
 J=\begin{pmatrix}1&-1\\1&1\end{pmatrix}.
\]
It has projective order four. Thus four divides \(|G_N^+|\).

For completeness, a finite-order element of \(\mathrm{PGL}_2(\mathbb Q)\)
can have order only 1, 2, 3, 4, or 6. Apart from the identity, its two complex
eigenvalues have ratio a primitive root of unity \(\zeta\), and
\[
 \frac{\operatorname{tr}(M)^2}{\det M}
 =2+\zeta+\zeta^{-1}.
\]
The right side is both rational and an algebraic integer, hence an integer.
It lies in [0,4]. The five possibilities for \(\zeta+\zeta^{-1}\)
are -2, -1, 0, 1, 2, giving exactly the stated orders. Finite order excludes
a nontrivial Jordan block when the eigenvalue ratio is one.

A generator of the cyclic group \(G_N^+\) therefore has order four.
The positive subgroup consists precisely of the four phase rotations by
multiples of a right angle. The reflection
\[
 q\longmapsto(\bar z_0/z_0)\bar q
\]
also preserves \(Q_N\), because it corresponds to z to \(\bar z\) on
the original integral circle. It is rational in half-angle coordinates,
represented by \(H\longmapsto\bar z_0\bar H\). Any other negative
element differs from it by an element of \(G_N^+\). Consequently
\[
 \boxed{G_N\cong D_4\text{ has order eight, and all its elements are
 conformal or anticonformal.}} \tag{1}
\]
This proof works even though the phase set need not have evenly spaced
points. Its exact quarter-turn symmetry and rationality are the essential
inputs.

## 2. Record norms force every nonincreasing-radius image to be a selfmap

Call N a strict record norm if
\[
 r_2(N)>\max_{1\le n<N}r_2(n).
\]
The set \(S_N\) is Gaussian primitive: otherwise dividing its common
Gaussian gcd would realize the same number of distinct integral points at
a smaller positive integer norm, contradicting the record property.
In particular its least squared realization radius is N.

Let a rational nonsingular projective map send its phases to a configuration
with least squared radius \(N'\le N\). Its image has exactly \(r_2(N)\)
distinct phases. An integral realization at N' has at most \(r_2(N')\)
points, so the record property forces \(N'=N\). Equality of cardinalities
then forces that realization to be the whole set \(S_N\).

The target anchor may differ from the source anchor. Their ratio is a
rational unit-circle number, so a common rational circle rotation identifies
the two normalized full phase sets. Such a rotation has a rational
half-angle representative. After this output normalization the original
map belongs to \(G_N\), and (1) shows that it is conformal or anticonformal.
Thus
\[
 \boxed{N'\le N\quad\Longrightarrow\quad M\text{ is conformal or
 anticonformal, for every full source at a strict record norm.}} \tag{2}
\]

Strict record cardinalities are unbounded. For example the Gaussian prime
factorization of \(5^k\) supplies the \(4(k+1)\) distinct representations
\(u(2+i)^j(2-i)^{k-j}\), with \(u\in\{1,-1,i,-i\}\) and
\(0\le j\le k\). Hence \(r_2(n)\) is unbounded, and so are its successive
strict record values. Formula (2) gives arbitrarily large configurations
with no nonconformal radius-nonincreasing map.

## 3. Limitation for the endpoint problem

Every full circle set contains a four-point quarter-turn orbit. Any connected
arc containing that orbit has angular width at least \(3\pi/2\). Thus
its endpoint constant is at least \((3\pi/2)N^{1/4}\), which diverges
along the record sequence. These examples do not meet the fixed endpoint
constant assumed in the proposed descent program.

The conclusion is that an existence theorem must use the short-arc
geometry in an essential way. Large cardinality, Gaussian primitivity,
and bounded-radius minimization alone cannot supply a nonconformal map.
An independent Luna audit confirmed the argument, including the distinction
between half-angle directions and ordinary radial lines and the required
target-anchor rotation.

The [exact symmetry checker](check_full_circle_record_symmetry.py)
reconstructs every possible cyclic or reversing permutation from three
labels and checks it on the full sets at norms 1, 5, 13, 25, 65, 85, and
325. Exactly eight maps survive in each case; their matrices satisfy the
conformal or anticonformal identities directly. It also verifies the
strict records through norm 325. These finite checks supplement the proof;
they do not establish the unbounded record statement.
