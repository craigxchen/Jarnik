# Direct common lifts: exact content cost and a sharp primitive four-point barrier

This note concerns common Gaussian multiplication of an actual source tuple;
it does not assume an admissible Möbius map exists. It gives an exact
denominator-clearing formula and a sharp power obstruction on the primitive
four-point family from [the Pell construction](four_point_bonus_counterexample.md).
It does not establish or refute any analogous assertion requiring eight or
more points, and does not prove the uniform endpoint bound.

Here a multiplier has **length** \(Q=|\mu|\), not norm \(Q^2\). A strip
width \(w\) means the diameter of the real coordinates of the multiplied
points. This is an unnormalized horizontal (radial) coordinate diameter. It
is different from an imaginary-coordinate strip normalized by the fourth
root of a squared radius. In particular the lower bound on \(w\) below is
not a theorem ruling out such a normalized strip or a theorem about Temur's
range of applicability.

## 1. Primitive tuples allow no hidden rational multiplier saving

Let \(z_0,\ldots,z_{m-1}\in\mathbf Z[i]\) have Gaussian gcd a unit. If
\(\mu\in\mathbf Q(i)\), then
\[
 \mu z_i\in\mathbf Z[i]\quad\hbox{for every }i
 \quad\Longleftrightarrow\quad \mu\in\mathbf Z[i].                 \tag{1}
\]
Indeed, Gaussian Bézout coefficients \(a_i\) satisfy \(\sum a_i z_i=1\),
so \(\mu=\sum a_i(\mu z_i)\) is integral. Moreover,
\[
 \gcd_i(\mu z_i)=\mu\quad\hbox{up to a Gaussian unit}.              \tag{2}
\]
Thus dividing the full common Gaussian content of a common integral lift
recovers the original tuple up to a unit.

There is a useful exact form for proposed phase lifts. Put
\(q_i=z_i/z_0\), fix \(0\ne\Gamma\in\mathbf Z[i]\), and seek the least
positive rational integer \(h\) such that every \(h\Gamma q_i\) is
Gaussian integral. Define
\[
 g=\gcd(z_0,\Gamma),\qquad A=z_0/g=x+iy,
 \qquad d=\gcd(|x|,|y|).
\]
Associates do not change the following quantities. Then
\[
 \boxed{\ h_{\min}=\frac{N(A)}d.\ }                              \tag{3}
\]
By (1), the condition is \(z_0\mid h\Gamma\), equivalently \(A\mid h\).
Writing \(A=d(x_1+iy_1)\) with coprime coordinates, the condition is
\(d^2(x_1^2+y_1^2)\mid hx\) and \(d^2(x_1^2+y_1^2)\mid hy\).
Their simultaneous least positive solution is
\(h=d(x_1^2+y_1^2)=N(A)/d\). This proof includes the prime \(2\), zero
coordinates, and nonprimitive individual points.

In the especially relevant case \(z_0=g_0\Gamma\), formula (3) becomes
\[
 h_{\min}=\frac{N(g_0)}{\gcd(\Re g_0,\Im g_0)},\qquad
 \frac{h_{\min}\Gamma}{z_0}
 =\frac{\overline{g_0}}{\gcd(\Re g_0,\Im g_0)}.                  \tag{4}
\]
The least denominator-cleared lift therefore costs precisely the length
of the ordinary primitive Gaussian part of the complementary factor.

If the original radius is \(R\), a common multiplier of length \(Q\)
gives radius \(S=QR\), squared radius \(Q^2R^2\), and multiplies arc
length by \(Q\). Its endpoint constant consequently changes by
\[
 C_{\rm lift}=\sqrt Q\,C_{\rm source}.                           \tag{5}
\]
These factors must be retained before any strip theorem is applied.

## 2. An actual primitive four-point family

Use the notation and proved primitivity of the
[four-point Pell family](four_point_bonus_counterexample.md):
\[
 U+V\sqrt5=(9+4\sqrt5)^n,\qquad n\equiv1\pmod{10},\qquad
 U^2-5V^2=1,
\]
\[
 P=-1+2i,\quad A=(U+V)+2Vi,\quad
 B=2V+(U-V)i,\quad C=(2U+8V)+(3U+V)i,
\]
\[
 (z_0,z_1,z_2,z_3)=
 (-PABC,-\bar P A\bar B\bar C,-\bar P\bar A B\bar C,
 -\bar P\bar A\bar B C).
\]
Write \(a=N(A),b=N(B),c=N(C)\). The original radius satisfies
\[
 R^2=5abc,\qquad R\asymp V^3.                                   \tag{6}
\]
The entire tuple is Gaussian primitive, each individual point has
coprime rational coordinates, and the points lie in an arc of length
\(L\asymp V\). In particular \(L/\sqrt R\asymp R^{-1/6}\to0\).
The upper arc bound is proved in the cited construction; the matching
lower bound follows from the first nonzero chord below. Its coordinate
expansion gives exactly
\[
\begin{aligned}
 z_1-z_0&=32V-i(16U+16V),\\
 z_2-z_0&=(-2U+2V)+4Vi,\\
 z_3-z_0&=(6U+2V)-i(4U+16V).
\end{aligned}                                                   \tag{7}
\]

## 3. Every common integral lift obeys a power tradeoff

For a nonzero Gaussian multiplier \(\mu=p-iq\), let
\[
 Q=\sqrt{p^2+q^2},\qquad
 w=\max_i\Re(\mu z_i)-\min_i\Re(\mu z_i).
\]
Here \(p,q\) are arbitrary integers; no primitivity of \(\mu\) is needed.
Then
\[
 \boxed{\quad w\ge \frac{16V}{5Q}-\frac{4Q}{V}.\quad}           \tag{8}
\]
To prove this, the first chord in (7) gives
\[
 w\ge16\left|(2p-q)V-qU\right|.
\]
Put \(a_0=2p-q\), \(b_0=q\). Since
\(a_0^2-5b_0^2\) is a nonzero integer,
\[
 |a_0-\sqrt5b_0|
 \ge\frac1{|a_0|+\sqrt5|b_0|}\ge\frac1{5Q}.
\]
Also
\[
 0<\frac UV-\sqrt5=rac1{V^2(U/V+\sqrt5)}<\frac1{4V^2}.
\]
The triangle inequality proves (8). In particular, if \(Q\le V/2\),
\[
 w\ge\frac{2V}{Q}.                                              \tag{9}
\]
Consequently, for any prescribed constant \(W_0>0\),
\[
 w\le W_0\quad\Longrightarrow\quad
 Q\ge\min\left(\frac12,\frac2{W_0}\right)V
 \asymp_{W_0}R^{1/3}.                                           \tag{10}
\]
By (1), (8)--(10) apply to **every rational Gaussian common multiplier
whose full output tuple is integral**. There is no rational denominator
loophole. More generally, \(Q\le V^{1-\eta}\), with fixed
\(0<\eta\le1\), forces \(w\gg V^\eta\asymp R^{\eta/3}\).

The exponent in (10) is attained. Take
\[
 \mu=\bar B=2V-i(U-V),\qquad Q=\sqrt b\asymp V.
\]
Using (7) and the Pell equation gives the exact relative real coordinates
\[
 \bigl(\Re(\mu(z_i-z_0))\bigr)_{i=0}^3=(0,-16,0,-4).             \tag{11}
\]
Thus \(w=16\). This is the integral normal to the chord \(z_2-z_0\).
For this family the least multiplier length capable of producing a
bounded horizontal coordinate diameter is therefore \(\Theta(R^{1/3})\).
The corresponding radius cost is
\[
 S=QR\asymp R^{4/3},\qquad S^2\asymp R^{8/3}.                   \tag{12}
\]
At this cost the endpoint constant in (5) is of constant order, since
the original constant has order \(R^{-1/6}\).

## 4. All six bounded pair residues incur the same clearing scale

Consider first the exceptionally near-real pair factor
\[
 \Gamma=PAC=X+i,\qquad z_0=-B\Gamma.
\]
Although its imaginary part is exactly \(1\), \(\Gamma\) alone does
not give a common integral lift of all four source phases. The block
\(B\) has coprime rational coordinates, so (4) gives exactly
\[
 h_{\min}=b\asymp V^2\asymp R^{2/3},\qquad
 \frac{h_{\min}\Gamma}{z_0}=-\bar B.                            \tag{13}
\]
The actual least denominator-cleared near-square lift has
\[
 S^2=b^2N(\Gamma)=(bX)^2+b^2=R^2b.                              \tag{14}
\]
Thus its **displayed** near-square gap is
\[
 b^2\asymp V^4\asymp S,                                        \tag{15}
\]
even though the uncleared factor had gap \(1\). This statement concerns
the displayed square \((bX)^2\), not an assertion about the distance to
every other integer square.

The same order of cost holds for each of the six near-real pair factors
in the construction. Their imaginary parts are respectively
\((-8,1,-1,-1,-3,2)\); each has norm \(\asymp V^4\), and its
complementary factor in either corresponding anchor has norm
\(\asymp V^2\). Every complementary factor has coprime rational
coordinates: its Gaussian prime factors come from disjoint rational
prime supports and use only one orientation at each prime. Formula (4)
therefore gives
\[
 h_{\min}\asymp V^2,\qquad Q\asymp V,\qquad
 S\asymp V^4,\qquad (h_{\min}\Im\Gamma)^2\asymp V^4.           \tag{16}
\]
Selecting any of these actual pair midpoints does not remove the power
loss. Selecting any other rational axis that produces a bounded
horizontal coordinate diameter is covered by the universal bound (10).
Exact alignment of an actual point with the real axis is still costlier:
because that point has coprime coordinates, the least nonzero integral
multiplier effecting this has length \(R\).

## 5. Scope of the obstruction

The earlier [projection route](fresh_geometric_projection_route.md)
explains the geometric strip criterion; the
[multipoint continuation](multipoint_continuation.md) already records
primitive chord normals. The additional statements here are the exact
common-content formula (3) and the sharp power tradeoff (8)--(12) on a
Gaussian primitive family, including all rational common integral lifts.

They show that bounded primitive pair residues and a vanishing endpoint
constant do not by themselves supply a cheap common bounded-width
horizontal lift, even for an actual primitive source configuration.
They do not exclude a different lift with a useful normalized imaginary
strip, a different arithmetic use of (14), or a conclusion that exploits
eight or more points simultaneously. The remaining direct-route question
is precisely whether larger cardinality forces an improvement over these
four-point multiplier and denominator costs.

All coordinate identities, the universal projection estimate, the
width-\(16\) multiplier, and the six pair-factor clearing costs are checked
in [the persistent exact checker](check_direct_primitive_axis_lift.py).
