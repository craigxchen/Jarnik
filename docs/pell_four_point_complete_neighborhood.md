# A complete endpoint neighborhood of the four-point Pell family

The known shrinking four-point family has no additional lattice point on
the same circle within chord distance \(\sqrt R/4\) of its designated
anchor, once its Pell parameter V is at least 16. This includes arbitrary
switches of individual Gaussian prime factors. It is stronger than the
earlier whole-block and reflection-orbit exclusions, but remains a theorem
about circles containing this specific four-point family. It does not
prove a bound for arbitrary endpoint clusters.

## 1. Statement and an integral coordinate frame

Take positive integers \(U,V\) with \(U^2-5V^2=1\), \(V\ge16\), and put
\[
 P=-1+2i,\quad A=U+V+2Vi,\quad B=2V+(U-V)i,
 \quad C=2U+8V+(3U+V)i.
\]
The four points are
\[
 (z_0,z_1,z_2,z_3)=
 (-PABC,-\bar P A\bar B\bar C,-\bar P\bar A B\bar C,
 -\bar P\bar A\bar B C).
\]
Set
\[
 a=\operatorname{Norm}(A),\quad b=\operatorname{Norm}(B),\quad
 c=\operatorname{Norm}(C)=16a-3b,\quad R^2=5abc.
\]
Then
\[
 \boxed{w\in\mathbb Z[i],\ |w|=R,\ |w-z_0|\le\sqrt R/4
 \quad\Longrightarrow\quad w\in\{z_0,z_1,z_2,z_3\}.} \tag{1}
\]
In particular an arc of length at most \(\sqrt R/4\) containing z0 has
at most four lattice points. The four designated points need not all belong
to this neighborhood at the smallest parameters; they do for sufficiently
large V because their chord lengths have order V and \(\sqrt R\) has
order \(V^{3/2}\).

Identify Gaussian integers with column vectors in \(\mathbb Z^2\), and define
\[
 T=\begin{pmatrix}2V&V-U\\-U-V&2V\end{pmatrix}.
\]
The Pell equation gives \(\det T=-1\). Direct expansion of the four
points yields
\[
 T^{-1}(z_i-z_0)=(0,0),(16,0),(0,2),(4,-6). \tag{2}
\]
Thus every lattice point w, regardless of its Gaussian factorization, has
integer coordinates \((s,t)=T^{-1}(w-z_0)\).

The Gram matrix and anchor covector are
\[
 G=T^TT=\begin{pmatrix}a&b-a\\b-a&b\end{pmatrix},\qquad
 T^Tz_0=(-8a,-b),\qquad ab-(a-b)^2=1. \tag{3}
\]
In particular \(\gcd(a,b)=1\). No coprimality assumptions on the individual
moving Gaussian blocks are needed for this argument.

## 2. Two fixed quadratic polynomials control every point on the circle

The equation \(|z_0+T(s,t)|^2=|z_0|^2\) is exactly
\[
 aF(s,t)+bH(s,t)=0, \tag{4}
\]
where
\[
 F(s,t)=s(s-2t-16),\qquad H(s,t)=t(2s+t-2).
\]
Since a and b are coprime, every integral solution satisfies
\[
 b\mid F(s,t),\qquad a\mid H(s,t). \tag{5}
\]
Both polynomials vanish at the four coordinate pairs in (2). Conversely,
their common zeros are exactly those four pairs: combine either
\(s=0\) or \(s=2t+16\) with either \(t=0\) or \(2s+t=2\).

Consequently it suffices to prove that a point in (1) has
\(|F(s,t)|<b\) and \(|H(s,t)|<a\). This retains exact concyclicity;
an approximate circle equation would not imply the divisibilities (5).

## 3. The endpoint neighborhood has small coordinates in the moving frame

For V at least four, the Pell equation implies
\[
 14V^2<a<15V^2,\qquad 5V^2<b<6V^2,\qquad c>206V^2,
 \qquad R>250V^3. \tag{6}
\]
For example \(2V<U\le9V/4\), and \(U>11V/5\); substituting these
in \(a=1+10V^2+2UV\), \(b=1+10V^2-2UV\) proves the first two bounds.
The last two follow from \(c=16a-3b\) and \(R^2=5abc\).

Let J be rotation by a right angle. From (3),
\[
 y_0=T^{-1}z_0=(b^2-9ab,\ 7ab-8a^2),\qquad
 v_0=T^{-1}Jz_0=(-b,8a).
\]
Thus
\[
 \|y_0\|_\infty<1800V^4,\qquad
 \|v_0\|_\infty<120V^2. \tag{7}
\]

Suppose \(|w|=|z_0|=R\), and write \(D=|w-z_0|\). Equal norms give
the exact radial component
\[
 w-z_0=-\frac{D^2}{2R^2}z_0+\beta Jz_0,
 \qquad |\beta|\le D/R.
\]
If \(D\le C_0\sqrt R\), applying \(T^{-1}\), (6), and (7) gives
\[
 \max(|s|,|t|)
 <\frac{18}{5}C_0^2 V+8C_0\sqrt V. \tag{8}
\]
For \(C_0=1/4\) and \(V\ge16\), the right side is at most
\[
 \frac9{40}V+2\sqrt V\le\frac{29}{40}V<V.
\]
It follows that
\[
 |F(s,t)|\le V(3V+16)\le4V^2<b,\qquad
 |H(s,t)|\le V(3V+2)<14V^2<a. \tag{9}
\]
Equations (5) and (9) force \(F=H=0\), proving (1).

## 4. What this resolves and what it leaves open

The [earlier fifth-point symmetry audit](four_point_pell_fifth_symmetry_obstruction.md)
excluded whole-block switches and points in the pair-reflection orbit, but
left switches of individual moving prime factors uncontrolled. Formula (1)
now covers every lattice point on the same circle in the stated
neighborhood, including all such switches. It also covers any fifth point
expressed with a different rational cotangent denominator; no fixed L=24
chart is assumed.

The proof does not transfer unchanged after multiplication by a nonunit
Gaussian integer. The transformed frame then has determinant of absolute
value greater than one, and a new lattice point can have nonintegral
coordinates in that frame. Clearing those denominators changes (5) and the
comparison in (9). Nor has an arbitrary large endpoint cluster been shown
to contain a quadruple of this fixed affine type. Those are separate
remaining issues; neither is supplied by this family theorem.

The [exact checker](check_pell_complete_neighborhood.py) verifies the moving
frame, the two fixed quadratic forms, and the displayed absolute constants
on 30 Pell instances. It also enumerates the entire integer circles for
the first four instances, containing 64, 256, 1,536, and 768 points. All
2,624 points obey the pencil equation and the neighborhood exclusion.
The larger circles include every partial prime-factor switch, not merely
the original whole blocks. These finite checks supplement the general
proof above.
