# Four common-unit points with bounded residues and vanishing endpoint constant

The four-point residue-growth target is false, including on primitive
circles, with four equal Gaussian units, and with endpoint constants
tending to zero. The following explicit infinite family has
\[
T=\prod_{i<j}t_{ij}=48,\qquad
V_1=\log5,\qquad V_2=W-\log5\longrightarrow\infty.
\tag{1}
\]
In particular no fixed \(\delta>0\) and \(B\) can make
\(2\log T\geq\delta V_2-B\) hold for every such four-point cluster.
This does not refute the analogous proposed eight-point target.

## 1. The exact family

Let
\[
U+V\sqrt5=(9+4\sqrt5)^n,\qquad n\equiv1\pmod{10}.
\tag{2}
\]
Thus \(U,V\) are positive integers,
\[
U^2-5V^2=1,\quad \gcd(U,V)=1,\quad U\equiv1\pmod4,\quad
V\equiv0\pmod4,\quad V\geq4.
\]
Define Gaussian integers
\[
\begin{aligned}
P&=-1+2i,\\
A&=(U+V)+2Vi,\\
B&=2V+(U-V)i,\\
C&=(2U+8V)+(3U+V)i.
\end{aligned}
\tag{3}
\]
The four points are
\[
\boxed{
\begin{aligned}
z_0&=-PABC,\\
z_1&=-\overline P\,A\overline B\,\overline C,\\
z_2&=-\overline P\,\overline A B\overline C,\\
z_3&=-\overline P\,\overline A\,\overline B C.
\end{aligned}}
\tag{4}
\]
All have the common norm
\[
N=R^2=5abc,\qquad
\begin{aligned}
a=N(A)&=1+10V^2+2UV,\\
b=N(B)&=1+10V^2-2UV,\\
c=N(C)&=13+130V^2+38UV.
\end{aligned}
\tag{5}
\]

For an explicit coordinate form, put
\[
K=(40UV^2+32U^2V)+i(60UV^2+4U^2V).
\]
Expansion of (4), using \(U^2-5V^2=1\), gives
\[
\begin{aligned}
z_0&=K+(U-17V)+i(8U+6V),\\
z_1&=K+(U+15V)-i(8U+10V),\\
z_2&=K-(U+15V)+i(8U+10V),\\
z_3&=K+(7U-15V)+i(4U-10V).
\end{aligned}
\tag{6}
\]
These formulas also exhibit their common large center and small
relative displacements.

## 2. Coprimality and the genuine common-unit representation

Each of \(A,B,C\) has relatively prime rational coordinates.
For \(A,B\), use \(\gcd(U,V)=1\), odd \(U\), and even \(V\).
For \(C\), the coordinate gcd divides \(22\); it is odd, and
divisibility by \(11\) would imply \(U\equiv-4V\pmod{11}\), contrary
to \(U^2-5V^2=1\pmod{11}\). All three norms are odd.
Consequently each block is coprime in \(\mathbf Z[i]\) to its
conjugate, and its rational prime divisors are all \(1\pmod4\).

The norms \(a,b,c\) are pairwise coprime. Indeed, \(a,b\) are
coprime to \(UV\), since each is congruent to \(-1\) modulo \(U\)
and to \(1\) modulo \(V\). Use
\[
a-b=4UV,\qquad c-13a=12UV,\qquad c-13b=64UV.
\tag{7}
\]
The only possible extra prime in the middle identity is \(3\);
it cannot divide the primitive sum of two squares \(a\).
The parity conditions exclude \(2\).

Modulo \(5\), the Pell recurrence implies
\[
UV\equiv n\pmod5.
\]
For example, the next product is congruent to \(UV+U^2=UV+1\).
Since \(n\equiv1\pmod5\), the three norms in (5) are respectively
\(3,4,1\pmod5\). Thus \(5,a,b,c\) are pairwise coprime.
In particular there is no overlap or prime-power adjustment at \(5\)
on the chosen progression.

There is a separate issue in passing from the block factorization (4)
to prime generators: arbitrary Gaussian blocks can contribute unit
factors, so a displayed common block unit alone would not suffice.
Here those factors can all be absorbed rigorously.

For odd \(n\), the Pell recurrence modulo \(3\) gives
\[
U\equiv0\pmod3,\qquad V\equiv\pm1\pmod3.
\]
Hence \(a,b,c\equiv2\pmod3\), so each is a nonsquare. Each block
therefore has a Gaussian prime occurring to an odd exponent \(e\).
Given a factorization \(H=i^j\prod\pi^{e_\pi}\), replace one such
prime generator by \(i^k\pi\), where \(ke\equiv j\pmod4\).
This makes the block factorization have unit \(1\). The norm supports
of \(A,B,C\) are disjoint, so these adjustments can be made
independently; choose \(P\) itself as the generator over \(5\).
Choose the conjugate generator over the conjugate prime throughout.

With these choices, every block in (3) is exactly a product of its
chosen prime generators without an extra unit. Consequently all four
points in (4) have the same prime-level Gaussian unit \(-1\).
This is stronger than merely having the same unit parity.

## 3. Primitive conductor and six exact residues

The four Gaussian integers have common gcd a unit. Every prime
factor of a block is allocated to both conjugate orientations somewhere
among the four rows, and no two different blocks share a norm prime.
Thus no common Gaussian normalization changes the family.

The prime over \(5\) has orientation counts \(1\) and \(3\).
Every prime-power layer from \(A,B,C\) has orientation counts \(2\)
and \(2\). Consequently
\[
W=\log N,\qquad V_1=\log5,\qquad
V_2=\log a+\log b+\log c=W-\log5.
\tag{8}
\]
This remains exact even when \(a,b,c\) have repeated prime factors.

For each pair write \(z_i=g_{ij}q_{ij}\) and
\(z_j=g_{ij}\overline{q_{ij}}\), including a shared sign in
\(g_{ij}\). The following choices are primitive, by the coprimality
already proved:

| Pair | \(q_{ij}\) | \(\operatorname{Im}q_{ij}\) |
|---|---|---:|
| \(0,1\) | \(PBC\) | \(-8\) |
| \(0,2\) | \(PAC\) | \(1\) |
| \(0,3\) | \(PAB\) | \(-1\) |
| \(1,2\) | \(A\overline B\) | \(-1\) |
| \(1,3\) | \(A\overline C\) | \(-3\) |
| \(2,3\) | \(B\overline C\) | \(2\) |

Each displayed imaginary part follows by expansion and
\(U^2-5V^2=1\). They are nonzero, so all four points are distinct.
Because \(q_{ij}\) and its conjugate are coprime, these are precisely
the primitive pair cofactors, up to harmless common signs. Therefore
\[
(t_{01},t_{02},t_{03},t_{12},t_{13},t_{23})
=(8,1,1,1,3,2),\qquad T=48.
\tag{9}
\]
The exact chord identity is
\[
|z_i-z_j|^2
=4t_{ij}^{\,2}R^2/N(q_{ij}).
\]
It uses the actual common units, not an artificial relaxation.

## 4. Endpoint constants tend to zero

As \(n\to\infty\) through the progression in (2),
\(U/V\to\sqrt5\). Formula (5) gives \(R\asymp V^3\).
Every displacement from \(K\) in (6) is \(O(V)\), so the four
points have diameter \(O(V)\).

For a direct uniform bound, \(U\leq3V\) implies that every displacement
from \(K\) has modulus at most \(40V\). Thus the diameter is at most
\(80V\). Also (5) gives
\[
a\geq10V^2,\qquad b\geq4V^2,\qquad c\geq130V^2,
\qquad R\geq\sqrt{26000}\,V^3.
\]
Relative to the argument of \(z_0\), every other point has smaller
angular distance at most \(2\arcsin(80V/(2R))\).
The bound \(\arcsin x\leq2x\) in this range places all four points
in an arc of length at most \(320V\). Hence they lie in an arc
of length \(C_n\sqrt R\), where
\[
C_n\leq
\frac{320}{26000^{1/4}}V^{-1/2}\longrightarrow0.
\tag{10}
\]
In particular every fixed positive endpoint constant, including
\(1/2\), contains these four-point examples for all sufficiently
large indices in the progression.

Combining (8)--(10),
\[
\frac{2\log T}{V_2}
=\frac{2\log48}{\log a+\log b+\log c}
\longrightarrow0.
\tag{11}
\]
Thus the proposed four-point residue bonus fails for every positive
coefficient, with no issue caused by a fixed large arc constant.

## 5. What the counterexample does and does not rule out

The exact four-point defect is \(D=V_1=\log5\).
At constant \(1/2\), the established inequality
\(D+2\log T\leq W+12\log(1/4)\) is entirely consistent with
\(W\to\infty\) and \(T=48\). A bounded unbalanced conductor can
therefore support arbitrarily much balanced conductor while keeping
all six primitive imaginary residues bounded.

The construction has only four points. It refutes attempts to prove
the eight-point endpoint theorem by imposing this stronger
four-point residue-growth inequality on every four-subcluster.
It neither constructs unbounded clusters nor refutes the eight-point
bonus or the uniform endpoint conjecture.

The common-unit condition here uses the explicitly permitted choice
of prime generators in Section 2. If prime generators were instead
required to be canonical primary associates, these four rows would
have units \(+1,+1,-1,-1\). Thus the example disproves a bound
required uniformly for common-unit models with arbitrary permissible
generators. The independent
[canonical-primary counterexample](canonical_four_point_counterexample.md)
closes the remaining unit convention issue: it has a different fixed
residue product, the same weights \(V_1=\log5\) and
\(V_2=W-\log5\), canonical common units, and endpoint constants
tending to zero. Both formulations of the four-point growth target
are therefore false.

A reflection construction must also check that its points are
distinct. For the canonical four-point family in the companion note,
the rows \(z_0=\overline F ABC\) and
\(z_3=FA\overline B\,\overline C\) satisfy
\[
\overline A z_3=\overline{\overline A z_0}.
\]
Thus multiplying that family by \(\overline A\) and adjoining all
conjugates produces at most six distinct points. Treating it as an
eight-row example would introduce duplicate rows and zero pair
residues; it supplies no eight-point counterexample.

## 6. A fixed affine relation survives the shrinking arc

This family also rules out a blanket lower bound for the height of an
integer affine relation among short-arc points. Write
`d_j = z_j - z_0` and identify Gaussian integers with their two real
coordinates. Formula (6) gives

\[
\begin{aligned}
d_1&=(32V,-16U-16V),\\
d_2&=(-2U+2V,4V),\\
d_3&=(6U+2V,-4U-16V).
\end{aligned}
\]

Using `det((x,y),(s,t)) = xt-ys` and `U^2-5V^2=1`, their three
anchored determinants are exactly

\[
\det(d_1,d_2)=-32,\qquad
\det(d_1,d_3)=96,\qquad
\det(d_2,d_3)=8.
\]

The two-dimensional determinant identity therefore gives
`8 d_1 - 96 d_2 - 32 d_3 = 0`, or the primitive affine relation

\[
\boxed{15z_0+z_1-12z_2-4z_3=0.}
\]

Its coefficients sum to zero and have maximum absolute value fifteen,
independently of the radius. The nonzero determinants show that the
four coordinate columns `(1, Re z_j, Im z_j)` have rank three, so this
relation generates their integer kernel. By Section 4, it persists
while `R` tends to infinity and the arc length divided by `sqrt(R)`
tends to zero. Consequently a proposed positive-power lower bound
for every nonzero affine-relation height fails already on four points;
even the existence of one bounded affine relation is not an endpoint
obstruction. This does not address a condition requiring more
independent relations, a growing number of rows, or the extracted full
cut profile. No uniform point-count conclusion follows from this audit.
