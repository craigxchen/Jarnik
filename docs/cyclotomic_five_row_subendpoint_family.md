# Five resonant Fibonacci rows in a sub-endpoint arc

This fixed one-class template produces five distinct lattice points on
one origin-centered circle whose five-point cluster has diameter
\(o(\sqrt R)\) as the common radius \(R\) grows.
It disproves a proposed four-row bound for arbitrary
mixed-prime cyclotomic layers. It does not give an unbounded point count.
The layer and height conventions are those of
[the resonant Fibonacci reduction](resonant_fibonacci_cyclotomic_reduction.md).

## Four odd orders and their divisor union

Let \(k\geq19\) be odd with \(k\equiv0\) or \(1\pmod3\), and put
\[
x=k-2,\quad y=k,\quad z=k+2,\quad w=k+4,
\qquad
n_1=xzw,\quad n_2=y^2w,\quad n_3=yz^2,\quad n_4=yzw,
\quad h=n_1.                                             \tag{1}
\]
The four numbers \(x,y,z,w\) are pairwise coprime. They are odd,
their differences are at most six, and the only possible nontrivial
pair gcd is \(\gcd(x,w)=3\), excluded by the congruence on \(k\).
Direct subtraction shows \(n_1<n_2,n_3<n_4\), and
\(n_4/h=y/x<3\). Hence every proper odd divisor of any \(n_j\)
lies below \(h\).

Let \(E=\bigcup_{j=1}^4\{e:e\mid n_j\}\), and set
\(U=\sum_{e\in E}\phi(e)\). Inclusion-exclusion, using
\(\sum_{e\mid n}\phi(e)=n\), yields
\[
U=\sum_{j=1}^4 n_j-(zw+yw+yz)
 =4k^3+15k^2-4k-24.                                    \tag{2}
\]
The only indices in \(E\) at least \(h\) are \(n_1,n_2,n_3,n_4\).
At \(k=21\), these are \(10925,11025,11109,12075\), and
\(U=43551\).

## Five orientations and their exact height

For \(e\in E\), put \(v_j(e)=\mathbf1_{e\mid n_j}\). Use
\[
v_1,\quad v_2,\quad v_3,\quad v_4,\quad
v_5=v_1+v_2+v_3.                                      \tag{3}
\]
Their base-layer coordinates are \(1,1,1,1,3\), so all
base-layer differences are even. For \(j=1,2,3,4\), the
Möbius coordinates satisfy
\[
\sum_{e:d\mid e}\mu(e/d)v_j(e)
 =\sum_{e:d\mid e\mid n_j}\mu(e/d)
 =\mathbf1_{d=n_j}.                                    \tag{4}
\]
The transform of \(v_5\) is the sum of those of
\(v_1,v_2,v_3\).
Thus every pair difference has zero coordinates at every odd
\(d<h\). Some pair, such as \(v_1,v_2\), first differs at \(h\).
All five vectors are distinct, and the minimum pair contact is
exactly \(h\), by [the odd-moment formula](resonant_fibonacci_cyclotomic_reduction.md#5-exact-odd-moment-constraints-and-their-limit).

Let \(W_e=\max_i v_i(e)-\min_i v_i(e)\) and
\(L=\sum_{e\in E}\phi(e)W_e\). For the first four rows the
width is one at every \(e\in E\setminus\{1\}\), and zero at
\(e=1\); their degree is \(U-1\). Adding \(v_5\) increases
the width at \(e\) by \((c_e-1)_+\), where
\(c_e=\sum_{j=1}^3\mathbf1_{e\mid n_j}\). This includes
\(e=1\), where the width rises from zero to two.
The identity \((c-1)_+=\binom c2-\binom c3\) for
\(c\in\{0,1,2,3\}\) gives
\[
\begin{aligned}
L&=U-1+\gcd(n_1,n_2)+\gcd(n_1,n_3)+\gcd(n_2,n_3)
       -\gcd(n_1,n_2,n_3)\\
 &=U-1+w+z+y-1
  =U+3k+4
  =4k^3+15k^2-k-20.                                 \tag{5}
\end{aligned}
\]
Consequently
\[
4h-L=k^2-15k-44>0\qquad(k\geq19).                   \tag{6}
\]
At \(k=21\), \(L=43618\) while \(4h=43700\).

## Equal-norm Gaussian realization

For odd \(d\), define the Gaussian Fibonacci factor
\(G_d=f_{(d+1)/2}+if_{(d-1)/2}\), where \(f_j\) is the
Fibonacci sequence. Choose \(d=4m+1\to\infty\), put
\(H_j=G_{n_jd}\), and let \(F=1+2i\). Define
\[
\begin{array}{ll}
Z_1=FH_1\overline H_2\overline H_3\overline H_4,
&Z_2=F\overline H_1H_2\overline H_3\overline H_4,\\
Z_3=F\overline H_1\overline H_2H_3\overline H_4,
&Z_4=F\overline H_1\overline H_2\overline H_3H_4,\\
Z_5=\overline F H_1H_2H_3\overline H_4.
\end{array}                                             \tag{7}
\]
All five have the same norm, since
\(N(F)=N(\overline F)=5\).
The factorization \(G_{n_jd}=\prod_{e\mid n_j}C_e\)
shows that their raw layer orientations are exactly (3).
The fifth row has two more oriented base factors than the
first four. Its prefactor ratio \(\overline F/F\) cancels the
fixed limiting-phase difference; this is the explicit
choice in [the converse construction](resonant_fibonacci_cyclotomic_reduction.md#6-converse-arbitrary-finite-layer-orientations-lift-to-a-template).
Therefore the five limiting directions coincide.

Let \(D_m\) be the Gaussian gcd of the \(Z_i\), and set
\(P_i=Z_i/D_m\). They have a common radius \(R_m\).
The layerwise common factor removes all growing common
content; the true gcd differs from it by bounded norm,
as proved in [the height reduction](resonant_fibonacci_cyclotomic_reduction.md#3-effective-primitive-radius-for-arbitrary-fixed-affine-rates).
It follows that
\[
\log R_m=\frac{dL\log\varphi}{2}+O(1),\qquad
R_m\longrightarrow\infty,                            \tag{8}
\]
where \(\varphi=(1+\sqrt5)/2\). Equation (4) and the
odd-moment expansion give every pair an angular difference
\(O(\varphi^{-dh})\), with some pair having nonzero leading
order \(h\). Thus
\[
\max_{i,j}|P_i-P_j|=O(R_m\varphi^{-dh}),\qquad
\frac{\max_{i,j}|P_i-P_j|}{\sqrt{R_m}}
 =O\!\left(\varphi^{-d(h-L/4)}\right)\longrightarrow0. \tag{9}
\]
The points are distinct for all sufficiently large \(m\).
For each fixed \(C>0\), five of them therefore lie in one
arc of length \(C\sqrt{R_m}\) for sufficiently large \(m\).
This proves a lower bound of five for any proposed uniform
point-count constant \(M(C)\); it does not settle whether
such a constant exists.

The four high indices also make the truncated Möbius matrix
on the divisor-closed set \(E\) have nullity four: its low
columns and rows form a triangular matrix with diagonal one.
Nullity alone does not give the fifth endpoint row. The even
base parity and the strict height budget (5)--(6) complete
the construction.

The dependency-free [exact checker](check_cyclotomic_five_row_subendpoint.py)
verifies the divisor identities, all odd moments through \(h\),
the Gaussian gcd, equal primitive norms, and a sufficient rational
inequality for an arc shorter than \(\sqrt R\) at the stated fixtures.
