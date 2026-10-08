# An exact three-row Gaussian construction with every nonempty cut

This note gives four actual Gaussian lattice points on one circle, with a full
seven-block incidence pattern among three nonanchor rows. Their six primitive
pair residues are bounded as the parameter grows. Their arc length is
\(O(R^{3/7})\), so the arc length divided by \(\sqrt R\) tends to zero.

The number of points is fixed at four. This is not a counterexample to the
uniform-count conjecture. It rules out an obstruction that would exclude even
three nonanchor rows with this balanced factor geometry.

## 1. Seven explicit Gaussian linear factors

For an integer parameter \(n\equiv2\pmod4\), put

\[
\begin{array}{c|c|c}
\text{block}&H(n)&\text{rows containing the block}\\ \hline
A&n-i&123\\
B&n-2-2i&12\\
C&3n-8+3i&13\\
D&3n-2+6i&1\\
E&3n+8+i&23\\
F&6n-11+16i&2\\
G&3n-16-i&3
\end{array}
\tag{1}
\]

Define

\[
h_1=ABCD,\qquad h_2=ABEF,\qquad h_3=ACEG.
\tag{2}
\]

Direct multiplication gives

\[
\begin{aligned}
h_1&=9n^4-48n^3+121n^2-152n+112+104i,\\
h_2&=18n^4-21n^3-8n^2+421n-26-442i,\\
h_3&=27n^4-144n^3-162n^2+944n-189-1088i.
\end{aligned}
\tag{3}
\]

Thus all three imaginary parts are fixed nonzero integers. The three pair
products after cancelling the common polynomial Gaussian factor also have
fixed imaginary parts:

\[
\begin{aligned}
EF\,\overline{CD}
 &=162n^4-405n^3+36n^2+3333n-6110-5850i,\\
EG\,\overline{BD}
 &=27n^4-144n^3-45n^2+632n-1840-1400i,\\
CG\,\overline{BF}
 &=54n^4-639n^3+2952n^2-7121n+7474-850i.
\end{aligned}
\tag{4}
\]

These are identities in \(\mathbb Z[i][n]\), not numerical approximations.

## 2. Actual points, actual radius, and an arc estimate

Let \(L=ABCDEFG\), let \(z_0=\overline L\), and for \(j=1,2,3\) set

\[
z_j=h_j\prod_{H\not\mid h_j}\overline H.
\tag{5}
\]

Every \(z_j\) lies in \(\mathbb Z[i]\), and all four have the same modulus

\[
R_n=|L(n)|=(486+o(1))n^7.
\tag{6}
\]

Also \(z_j/z_0=h_j/\overline{h_j}\). From (3),

\[
\arg(z_j/z_0)=2\arg h_j=O(n^{-4}).
\tag{7}
\]

Taking the branch near zero is valid for all sufficiently large positive
\(n\), because each real part in (3) has positive leading coefficient. All
four points therefore lie in one arc of angular width \(O(n^{-4})\), hence
of length

\[
O(n^3)=O(R_n^{3/7}),\qquad
\frac{\text{arc length}}{\sqrt{R_n}}=O(n^{-1/2})\longrightarrow0.
\tag{8}
\]

Every pair is distinct: (3) and (4) show that its half-angle numerator has
nonzero imaginary part for every integer parameter. For example,

\[
z_j-z_0=2i\operatorname{Im}(h_j)
                \prod_{H\not\mid h_j}\overline H
\tag{9}
\]

also directly gives the \(O(n^3)\) chord estimate.

## 3. Specialization, primitivity, and the exact scale of the blocks

The roots of the seven factors are

\[
i,\ 2+2i,\ \tfrac83-i,\ \tfrac23-2i,\
-\tfrac83-\tfrac13i,\ \tfrac{11}6-\tfrac83i,\
\tfrac{16}3+\tfrac13i.
\tag{10}
\]

They are distinct, none is real, and no two are complex conjugates. Thus any
two different factors, or a factor and the conjugate of any factor, have
nonzero constant resultants over \(\mathbb Q(i)\). Clearing the fixed
denominators gives a nonzero Gaussian integer which is divisible by every
common divisor of their values at an integer \(n\).

Consequently all unintended overlaps on specialization have uniformly
bounded Gaussian valuations and modulus. One may remove these bounded
overlaps into correction factors to obtain seven mutually coprime,
conjugate-coprime Gaussian blocks \(H_S'\), with

\[
\log N(H_S')=2\log n+O(1)
\quad(S\ne\varnothing),\qquad
\log|K_j|=O(1).
\tag{11}
\]

For completeness, at a Gaussian prime occurring in several evaluated blocks,
retain the largest valuation and remove the others. The removed valuations
are bounded by the fixed pair resultants. If both conjugate orientations
occur, retain the larger orientation and remove the bounded smaller one;
the cross-conjugate resultants bound this loss. Inert or ramified common
contents are bounded by the self-conjugate resultants. There are only
finitely many exceptional primes, all determined by the fixed factors.

The four polynomials in (5) have no common nonconstant Gaussian polynomial
factor: for each of the seven cuts, some row flips its orientation relative
to the anchor. Their common Gaussian divisor after specialization is
therefore bounded by a fixed nonzero Gaussian integer, by a polynomial
Bezout identity with cleared denominators. Dividing out their actual common
divisor gives a primitive tuple with intrinsic radius

\[
R_n^{\mathrm{prim}}\asymp n^7.
\tag{12}
\]

Its arc still satisfies (8) with this primitive radius. Equivalently, the
least Gaussian realization radius of these four rational directions has
order \(n^7\): their Gaussian lcm contains all seven polynomial blocks,
up to the same bounded resultant factors. Any alternative realization
differs by a common rational Gaussian multiplier, and a primitive integral
tuple permits only an integral multiplier.

The congruence \(n\equiv2\pmod4\) ensures that
\(A,B/2,C,D/2,E,F,G\) are integral and all have odd norms. Thus, after
removing ordinary integer coordinate contents from (3) and (4), the two
coordinates have opposite parities. The resulting Gaussian numerators are
conjugate-coprime, with bounded nonzero imaginary parts. In particular no
division by \(1+i\), and hence no unit rotation of a near-real numerator,
is needed. Such a rotation would not preserve the assertion that its
imaginary part is small. Additional pair-gcd corrections on specialization
are bounded by the same fixed resultants, so all six actual primitive pair
residues stay bounded.

There are seven nonempty blocks, not eight. Relative to the intrinsic
\(W=2\log R_n^{\mathrm{prim}}\sim14\log n\), their individual weights are
\(W/7+O(1)\), and the empty block has zero weight. If an inherited full
eight-block profile is desired, multiply all four points by the conjugate
of an independent linear Gaussian block \(H_\varnothing(n)\). The new
inherited radius is of order \(n^8\), each of the eight block weights is
\(2\log n+O(1)\), and the arc has order \(\sqrt R\). That empty factor is a
common Gaussian divisor and disappears on primitive reduction. One must
not identify the inherited and intrinsic values of \(W\).

## 4. How the rational compatibility was found

Use monic linear factors with roots

\[
a=i,\quad b=3s+2i,\quad c=4s-i,\quad d=s-2i,
\qquad s\in\mathbb Q.
\tag{13}
\]

Their product \(f_1\) has real nonconstant coefficients and constant
imaginary part \(12s(1+s^2)\). Put

\[
f_2=f_1-\lambda N((n-a)(n-b)),\qquad
f_3=f_1-\mu N((n-a)(n-c)).
\tag{14}
\]

For real \(\lambda,\mu\), all imaginary parts remain unchanged; the first
pair shares the factors with roots \(a,b\), and the first and third share
the roots \(a,c\). A new common root \(z=x+iy\) of the two remaining
quadratics must make

\[
\lambda=\frac{(z-c)(z-d)}{(z-\bar a)(z-\bar b)},\qquad
\mu=\frac{(z-b)(z-d)}{(z-\bar a)(z-\bar c)}
\tag{15}
\]

real. Away from the zero and pole cases, the quotient of these two
expressions first forces the rational circle

\[
x^2+y^2+(3/s-7s)x+12s^2-13=0.
\tag{16}
\]

Parametrize this circle by the line \(y+1=t(x-4s)\). Apart from its base
point, this gives

\[
x=\frac{3s+4st^2+2t-3/s}{1+t^2},\qquad
y=\frac{t^2-(s+3/s)t-1}{1+t^2}.
\tag{17}
\]

After the zero/pole factors are removed, the remaining real-compatibility
condition is

\[
2s(s^2+2)t^2+(s^2-3)t+s(s^2-1)=0.
\tag{18}
\]

Its discriminant is

\[
(s^2+1)^2(9-8s^2).
\tag{19}
\]

Thus the conic \(r^2+8s^2=9\) supplies rational parameters. Take
\(s=2/3\), \(r=7/3\), and \(t=-1/8\). Then

\[
z=-\tfrac83-\tfrac13i,\qquad
\lambda=\tfrac{25}{17},\qquad
\mu=\tfrac{175}{136}.
\tag{20}
\]

The other two roots are \(11/6-8i/3\) and \(16/3+i/3\). Clearing the
linear-factor denominators gives exactly (1)–(4). In particular the
construction does not depend on approximate root finding or an uncharged
algebraic coefficient field.

The reduction (18) was also checked as an exact bivariate polynomial
identity. If \(Z=s(1+t^2)(x+iy)\) and \(D_0=s(1+t^2)\), substitution of
(17) gives

\[
\begin{aligned}
&\operatorname{Im}\bigl[(Z-cD_0)(Z-dD_0)
       (\overline Z-aD_0)(\overline Z-bD_0)\bigr]\\
&\quad=6s^2(1+t^2)^2(st-1)(2st-s^2-3)\\
&\qquad\qquad\cdot
 \bigl(2s(s^2+2)t^2+(s^2-3)t+s(s^2-1)\bigr).
\end{aligned}
\tag{21}
\]

All 29 nonzero coefficients and the discriminant factorization (19) were
checked with exact integer polynomial arithmetic.

As a separate normalization check, the first 1,000 positive parameters
\(n\equiv2\pmod4\) were evaluated with exact Gaussian integer arithmetic.
All 6,000 pair numerators were conjugate-coprime after primitive reduction,
with nonzero residues. The largest observed magnitudes, in the order
\(01,02,03,12,13,23\), were \(26,221,1088,2925,350,425\). These finite
checks supplement the polynomial and resultant proofs above.

## 5. Consequence and remaining question

All one-row small-residue equations, all three cross-row overlap equations,
and the actual least-radius accounting coexist in this example. Merely
adding the third row to the one-row cyclotomic tests does not force a
positive block-scale correction or residue height.

The further question is whether another family can continue this compatibility
to four or more nonanchor rows while retaining every required cut and a degree
budget comparable with half of the total block degree. For the specific triple
above, a subsequent [exact exhaustive calculation](quartic_extension_maximality.md)
excludes any fourth constant-imaginary polynomial anchor, of any degree.
That family-specific obstruction does not exclude other polynomial systems
or isolated integer extensions.
