# Independent audit: the Pell four-point cluster is locally maximal

## Result

Let
\[
 U+V\sqrt5=(9+4\sqrt5)^n,\qquad n\equiv1\pmod {10},
 \qquad U^2-5V^2=1,
\]
and let \(z_0,z_1,z_2,z_3\) be the primitive Pell circle points from
[the four-point construction](four_point_bonus_counterexample.md).  Write
\[
 |z_i|=R.
\]
For every member with \(V\ge16\),
\[
 \boxed{\quad
 w\in\mathbb Z[i],\quad |w|=R,\quad
 |w-z_0|\le\frac14\sqrt R
 \ \Longrightarrow\
 w\in\{z_0,z_1,z_2,z_3\}.
 \quad}                                                   \tag{1}
\]
Thus the Euclidean ball of radius \(\sqrt R/4\) about \(z_0\) contains at
most the four displayed lattice points from that circle.  Some of the
four points need not lie in the ball for a small family member; assertion
(1) only classifies the points that do.

Every circular arc through \(z_0\) of length at most \(\sqrt R/4\) is
contained in this Euclidean ball in the needed pairwise sense, since chord
distance is at most arc distance.  Hence this particular Pell cluster
cannot be enlarged to five points in such an arc on the same primitive
circle.

The theorem is specific to this four-point family and this circle.  It does
not give a general five-point bound.  Multiplication by a nonunit Gaussian
integer changes the primitive lattice coordinates and is not covered by
the unimodular-coordinate argument below.

## 1. A unimodular coordinate system

Use
\[
\begin{aligned}
 P&=-1+2i,\\
 A&=(U+V)+2Vi,\\
 B&=2V+(U-V)i,\\
 C&=(2U+8V)+(3U+V)i
\end{aligned}
\]
and
\[
 (z_0,z_1,z_2,z_3)=
 (-PABC,-\bar P A\bar B\bar C,-\bar P\bar A B\bar C,
  -\bar P\bar A\bar B C).
\]
Set
\[
 a=N(A)=U^2+2UV+5V^2,\qquad
 b=N(B)=U^2-2UV+5V^2.
\]
The integer matrix
\[
 T=\begin{pmatrix}
 2V&V-U\\
 -U-V&2V
 \end{pmatrix}                                          \tag{2}
\]
has
\[
 \det T=5V^2-U^2=-1.                                    \tag{3}
\]
The four displacements have the exact coordinates
\[
 z_i-z_0=T(s_i,t_i),\qquad
 (s_i,t_i)\in\{(0,0),(16,0),(0,2),(4,-6)\}.             \tag{4}
\]
Indeed (4) expands to
\[
\begin{aligned}
z_1-z_0&=(32V,-16U-16V),\\
z_2-z_0&=(-2U+2V,4V),\\
z_3-z_0&=(6U+2V,-4U-16V).
\end{aligned}
\]

The Gram matrix and its determinant are
\[
 G=T^TT=
 \begin{pmatrix}a&b-a\\b-a&b\end{pmatrix},\qquad
 ab-(a-b)^2=1.                                          \tag{5}
\]
In particular \(\gcd(a,b)=1\).  Direct multiplication gives
\[
 T^Tz_0=(-8a,-b).                                       \tag{6}
\]

## 2. The exact integral conic

Let \(w\in\mathbb Z[i]\) lie on the same circle, so \(|w|=R=|z_0|\).
Because \(T\) is unimodular, there is a unique
\((s,t)\in\mathbb Z^2\) with
\[
 w-z_0=T(s,t).
\]
Expanding \(|z_0+T(s,t)|^2=|z_0|^2\) with (5)--(6) gives
\[
 \boxed{\quad
 aF(s,t)+bH(s,t)=0,\quad}                               \tag{7}
\]
where
\[
 F(s,t)=s(s-2t-16),\qquad
 H(s,t)=t(2s+t-2).                                      \tag{8}
\]
Since \((a,b)=1\),
\[
 b\mid F(s,t),\qquad a\mid H(s,t).                      \tag{9}
\]

The simultaneous zero set of the two quadratic factors is exactly
\[
 F=H=0
 \quad\Longleftrightarrow\quad
 (s,t)\in\{(0,0),(16,0),(0,2),(4,-6)\}.                 \tag{10}
\]
These are precisely the four old points in (4).

## 3. Quantitative control near \(z_0\)

The third block norm is
\[
 c=N(C)=16a-3b.
\]
For \(V\ge4\), the Pell equation gives
\[
 14V^2<a<15V^2,\qquad
 5V^2<b<6V^2,\qquad
 c>206V^2.                                              \tag{11}
\]
Since
\[
 R^2=5abc,
\]
these inequalities imply
\[
 R>250V^3.                                              \tag{12}
\]

Put
\[
 y_0=T^{-1}z_0,\qquad v_0=T^{-1}Jz_0,
\]
where \(J(x,y)=(-y,x)\).  From (5)--(6),
\[
\begin{aligned}
y_0&=\bigl(-(9ab-b^2),-(8a^2-7ab)\bigr),\\
v_0&=(-b,8a).
\end{aligned}                                           \tag{13}
\]
Equations (11) give the safe coordinate bounds
\[
 \|y_0\|_\infty<1800V^4,\qquad
 \|v_0\|_\infty<120V^2.                                 \tag{14}
\]

Let \(D=|w-z_0|\).  On the circle one can write, without choosing an
angular branch,
\[
 w-z_0=-\frac{D^2}{2R^2}z_0+\beta Jz_0,\qquad
 |\beta|\le\frac DR.                                    \tag{15}
\]
Applying \(T^{-1}\), (12)--(14) show that
\[
 M:=\max(|s|,|t|)
 \le 3.6\,C^2V+8C\sqrt V
 \qquad\text{whenever }D\le C\sqrt R.                   \tag{16}
\]
For \(C=1/4\) and \(V\ge16\),
\[
 M\le\frac9{40}V+2\sqrt V
 \le\frac{29}{40}V<V.                                   \tag{17}
\]

It follows from (8), (11), and \(V\ge16\) that
\[
\begin{aligned}
 |F(s,t)|
 &\le M(3M+16)<V(3V+16)\le4V^2<b,\\
 |H(s,t)|
 &\le M(3M+2)<V(3V+2)<14V^2<a.
\end{aligned}                                           \tag{18}
\]
The divisibilities (9) now force \(F=H=0\), and (10) proves (1).

## 4. Scope

This argument uses the extra integral equation for a point on the same
circle, together with a unimodular coordinate system adapted to the four
known displacements.  It closes the possible fifth-point extension of this
Pell family inside one fixed endpoint-scale neighborhood.

It does not constrain a fifth point obtained after changing the primitive
circle, applying a nonunit common Gaussian scaling, or starting from an
unrelated endpoint configuration.  In particular it supplies no general
bound for five or more lattice points and no improvement to the general
point-count growth rate.

The identities and all safe constants are verified in
[the independent exact checker](check_pell_four_point_local_ball_independent.py).
