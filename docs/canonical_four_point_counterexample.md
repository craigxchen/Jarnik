# A four-point counterexample with canonical primary units

The four-point balanced-conductor residue-growth target fails even
when every Gaussian prime generator is required to be its primary
associate. The following primitive infinite family has four points
with the same canonical unit, fixed primitive pair residues, and arc
length divided by the square root of the radius tending to zero.
It strengthens the unit audit in
[the smaller-residue counterexample](four_point_bonus_counterexample.md).
It does not disprove an eight-point bonus or the uniform endpoint
conjecture.

## 1. The family

Put
\[
\lambda=9+4\sqrt5,\qquad
x_r+y_r\sqrt5=\lambda^r,\qquad
F=1+2i,\qquad H_r=x_r+Fy_r.
\tag{1}
\]
For every nonnegative integer \(r\),
\[
x_r^2-5y_r^2=1,\qquad \gcd(x_r,y_r)=1,
\qquad x_r\equiv1\pmod4,\quad y_r\equiv0\pmod4.
\tag{2}
\]
The congruences follow from the recurrence
\[
x_{r+1}=9x_r+20y_r,\qquad y_{r+1}=4x_r+9y_r.
\]
Choose \(n\equiv1\pmod{10}\), \(n\geq1\), and set
\[
a=n,\quad b=n+2,\quad c=n+4,\qquad
A=H_a,\quad B=H_b,\quad C=H_c.
\]
Define
\[
\begin{aligned}
z_0&=\overline F\,ABC,\\
z_1&=F\overline A\,\overline B\,C,\\
z_2&=F\overline A\,B\overline C,\\
z_3&=FA\overline B\,\overline C.
\end{aligned}
\tag{3}
\]
Their common radius is
\[
R=\sqrt5\,|A|\,|B|\,|C|,\qquad N=R^2.
\tag{4}
\]

## 2. Primitive and disjoint block supports

Each \(H_r=(x_r+y_r)+2iy_r\) has relatively prime rational
coordinates: a common divisor must divide \(2\), and its real
coordinate is odd. Its norm is odd. Consequently \(H_r\) is
coprime to \(\overline H_r\); no prime \(2\), or rational prime
congruent to \(3\pmod4\), divides its norm.
Also \(H_r\) and the rational integer \(y_r\) are coprime in
\(\mathbf Z[i]\), by \(\gcd(x_r,y_r)=1\).

Here is a direct check that the three different blocks have disjoint
rational norm supports. Pell multiplication gives
\[
\begin{aligned}
H_{r+d}
 &=x_dH_r+y_d(5y_r+Fx_r),\\
\overline H_{r+d}
 &=x_d\overline H_r+y_d(5y_r+\overline F x_r).
\end{aligned}
\]
Reducing modulo \(H_r\), and using \(F\overline F=5\), yields
\[
H_{r+d}\equiv4(2-i)y_dy_r\pmod{H_r},\qquad
\overline H_{r+d}\equiv-4ix_dy_r\pmod{H_r}.
\tag{5}
\]
For the only differences needed here,
\[
\begin{array}{c|cc}
d&x_d&y_d\\ \hline
2&161=7\cdot23&72=2^3\,3^2\\
4&51841=47\cdot1103&23184=2^4\,3^2\,7\cdot23 .
\end{array}
\tag{6}
\]
All displayed odd prime factors are \(3\pmod4\); \(1103\) is
prime (trial division up to \(31\) suffices). Thus a common
Gaussian prime of \(H_r,H_{r+d}\), or of
\(H_r,\overline H_{r+d}\), could only lie over \(2\) or be
associated to \(2-i=-iF\). The first is already excluded, and
\(F\) never divides \(H_r\): modulo \(F\), the latter equals
\(x_r\), while \(x_r^2\equiv1\pmod5\).
This proves all required pairwise and cross-conjugate coprimality
among \(A,B,C\).

To check the remaining factor over \(5\), the recurrence modulo \(5\)
gives
\[
x_r\equiv(-1)^r,\qquad y_r\equiv(-1)^r r\pmod5.
\]
Modulo \(\overline F=1-2i\), one has \(F\equiv2\), so
\[
\overline F\mid H_r
\quad\Longleftrightarrow\quad r\equiv2\pmod5.
\tag{7}
\]
The three indices \(a,b,c\) are respectively \(1,3,0\pmod5\),
so this never occurs. Hence \(F,A,B,C\) have pairwise disjoint
rational norm supports and are individually conjugate-coprime.

Every block in (3) occurs in both conjugate orientations among
the four rows. Therefore the common Gaussian gcd of the four
points is a unit: the whole circle configuration is primitive.

## 3. The canonical units really agree

Use the convention that the primary associate of an odd Gaussian
integer is the associate congruent to \(1\) modulo
\((1+i)^3\). This convention selects one associate of each odd
Gaussian prime, and conjugation preserves being primary.

Equation (2) gives \(H_r\equiv1\pmod4\), so every \(H_r\) and
its conjugate are primary. Their factorizations into primary prime
generators have unit \(1\). On the other hand,
\[
F\equiv\overline F\equiv-1\pmod{(1+i)^3}.
\]
Thus \(F\) and \(\overline F\) both have canonical unit \(-1\).
All four rows in (3) consequently have the same canonical unit
\(-1\), with no choice of noncanonical associates.

Equivalently, all four displayed points are congruent to
\(1+2i\pmod4\), which also directly verifies equality of their
classes modulo \((1+i)^3\).

## 4. Exact fixed primitive residues

For \(s>r\), Pell subtraction gives
\[
x_rx_s-5y_ry_s=x_{s-r},\qquad
x_ry_s-x_sy_r=y_{s-r}.
\]
Expanding (1) therefore gives the two exact identities
\[
\operatorname{Im}(\overline F H_rH_s)=-2x_{s-r},\qquad
\operatorname{Im}(\overline H_rH_s)=2y_{s-r}.
\tag{8}
\]
For each pair in (3), factor \(z_i=gq\), \(z_j=g\overline q\).
The primitive cofactors and their absolute imaginary parts are:

| Pair | Cofactor \(q\) | \(t_{ij}=|\operatorname{Im}q|\) |
|---|---|---:|
| \(0,1\) | \(\overline F AB\) | \(2x_2=322\) |
| \(0,2\) | \(\overline F AC\) | \(2x_4=103682\) |
| \(0,3\) | \(\overline F BC\) | \(2x_2=322\) |
| \(1,2\) | \(\overline B C\) | \(2y_2=144\) |
| \(1,3\) | \(\overline A C\) | \(2y_4=46368\) |
| \(2,3\) | \(\overline A B\) | \(2y_2=144\) |

The coprimality proved above ensures \(\gcd(q,\overline q)=1\),
so these are the primitive pair cofactors, up to an irrelevant sign.
The use of canonical primary generators also changes a displayed
cofactor by at most a sign: cofactors involving \(F\) have unit
\(-1\), and the other cofactors have unit \(1\).
Thus the same absolute imaginary parts are the residues for the
canonical common-unit model. Every residue is nonzero, so the points
are distinct, and
\[
T=\prod_{i<j}t_{ij}
=2^6x_2^2x_4y_2^2y_4
=10336141769048653824
\tag{9}
\]
is independent of \(n\).

The prime over \(5\) is the only block with a \(1\)-versus-\(3\)
orientation cut. Every inherited prime-power layer from \(A,B,C\)
has a \(2\)-versus-\(2\) cut. Writing \(W=\log N\), the private and
balanced layer weights are exactly
\[
V_1=\log5,\qquad
V_2=\log N(A)+\log N(B)+\log N(C)=W-\log5.
\tag{10}
\]
Repeated primes inside an individual block cause no change to this
layer count.

## 5. The endpoint constants tend to zero

For \(r\geq1\), the positive Pell coordinates give
\[
\frac{\lambda^r}{2}\leq x_r\leq |H_r|
\leq x_r+\sqrt5\,y_r=\lambda^r.
\]
Consequently
\[
\frac{\sqrt5}{8}\lambda^{3n+6}
\leq R\leq\sqrt5\,\lambda^{3n+6}.
\tag{11}
\]
The three differences from \(z_0\) have exact lengths
\[
|z_0-z_1|=4x_2|C|,\quad
|z_0-z_2|=4x_4|B|,\quad
|z_0-z_3|=4x_2|A|.
\]
Each is at most \(L_n=4x_4\lambda^{n+4}\), with
\(L_n/R\to0\). For all sufficiently large \(n\), each point's
smaller angular distance from \(z_0\) is at most
\[
2\arcsin(L_n/(2R))\leq2L_n/R.
\]
They therefore lie in an arc of length at most \(4L_n\). By (11),
\[
\frac{\text{arc length}}{\sqrt R}
\leq
16x_4\sqrt{\frac8{\sqrt5}}\,
\lambda^{1-n/2}
\longrightarrow0.
\tag{12}
\]
Thus every fixed positive endpoint constant, including \(1/2\),
contains these four-point examples for all sufficiently large
\(n\equiv1\pmod{10}\).

Equations (9)--(12) disprove, for every \(\delta>0\), a uniform
inequality
\[
2\log T\geq\delta V_2-O(1)
\]
on primitive four-point endpoint clusters, even with canonical
primary common units.

## 6. The exact-balance limitation

The exact-balanced four-row rational-ray obstruction remains valid.
Here the fixed private factor has norm \(5\), so exact balance is
absent despite \(V_1/W\to0\). Indeed
\[
\lambda^{-r}H_r\longrightarrow
\frac{1+F/\sqrt5}{2}
\]
has slope
\[
\frac{2}{\sqrt5+1}=\frac{\sqrt5-1}{2},
\]
a quadratic irrational. Identity (8) says that the imaginary
remainders of products along this ray stay fixed, and
\(\overline F(1+F/\sqrt5)^2\) is real. This is the precise
mechanism by which one fixed private Gaussian factor supports
arbitrarily much balanced conductor and bounded residues.
It excludes a stability claim based only on the private conductor
being bounded, or being a vanishing proportion of total conductor.
It does not exclude a genuinely eight-row mixed-layer argument.

The derivations above use exact integer and Gaussian identities.
Numerical spot checks at \(n=1,11,21\) additionally confirmed the
norm equality, all six residues, the pairwise norm gcds, and the
common residue class modulo \(4\). No Lean formalization of this
counterexample is claimed.
