# A high-core exact graph can preserve an actual four-point endpoint cluster

This is a bounded obstruction to a proposed archimedean graph-rigidity
argument. Every four rational circle phases have an exact rational
projective involution interchanging opposite cyclic labels. On a shrinking
four-point endpoint cluster its invariant coefficient core is unbounded,
yet its target is precisely the source and its least radius is unchanged.
Thus a theorem forcing growth for all high-height cores needs more than
four points and more than the exact four-point interpolation constraints.
This does not obstruct a theorem with a sufficiently large number of points.

The source-capped and actual-anchor estimates in
[mobius_source_capped_content.md](mobius_source_capped_content.md) and
[mobius_actual_anchor_normalization.md](mobius_actual_anchor_normalization.md)
do not contradict this example. Its determinant need not be one, and its
large coefficient-core ratio cannot be removed by choosing another actual
source anchor.

## 1. Exact rational involution of any four points

After a common rational circle rotation, write four ordered half-angle
parameters as
\[
 t_0<t_1<t_2<t_3,\qquad t_i\in\mathbb Q.
\]
Put
\[
 c=t_0+t_2-t_1-t_3<0,\qquad
 a=t_0t_2-t_1t_3,\qquad
 b=ct_0t_2-a(t_0+t_2).
\]
Then
\[
 F(t)=\frac{at+b}{ct-a}
\]
interchanges \(t_0,t_2\) and interchanges \(t_1,t_3\). Indeed the
condition \(F(x)=y\) is
\(a(x+y)+b-cxy=0\), and the definitions impose this for both pairs.
The trace-zero matrix squares to a scalar, so each interchange works in
both directions.

To see that this map is nonsingular and orientation preserving, set
\(h=a/c\). Direct subtraction gives
\[
 h-t_1=\frac{(t_0-t_1)(t_2-t_1)}c>0,\qquad
 h-t_2=\frac{(t_2-t_1)(t_3-t_2)}c<0.
\]
Thus \(t_1<h<t_2\). Put
\(s^2=(h-t_0)(t_2-h)>0\). The same map becomes
\[
 F(t)=h-\frac{s^2}{t-h},\qquad
 A=\begin{pmatrix}h&-h^2-s^2\\1&-h\end{pmatrix},
 \qquad\det A=s^2>0. \tag{1}
\]
Clearing denominators and removing the common integer content produces a
primitive integral representative. The pole \(h\) lies between the middle
points, but no source point is a pole. The image as a finite phase set is
exactly the original four-point set.

If the phase arc has width less than \(\pi\), this map is nonconformal.
A positive-determinant conformal map is a rotation on phases. A nontrivial
involution of that kind is the antipodal rotation, which cannot preserve a
set contained in such an arc.

## 2. The invariant core must be large

Let \(\varepsilon=t_3-t_0\). Since
\(s^2=(h-t_0)(t_2-h)\),
\[
 s^2\le(t_2-t_0)^2/4\le\varepsilon^2/4.
\]
For (1), the scale-invariant trace-to-determinant ratio therefore satisfies
\[
 \frac{T}{\delta}
 =\frac{1+2h^2+(h^2+s^2)^2}{s^2}
 \ge\frac4{\varepsilon^2}. \tag{2}
\]
Here \(T=\operatorname{tr}(A^TA)\); the same ratio holds after clearing
rational denominators. Exchanging the two homogeneous coordinates to use
the other half-angle convention is orthogonal and does not affect it.

For the complex-linear/antilinear coefficients \(\alpha,\beta\),
\[
 r=\frac{N\beta}{N\alpha}
 =\frac{T-2\delta}{T+2\delta}\in(0,1).
\]
Write this invariant ratio in lowest terms as \(a_c/b_c\). Then
\[
 b_c\ge\frac1{1-r}
 =\frac{T+2\delta}{4\delta}
 \ge\varepsilon^{-2}+\frac12. \tag{3}
\]
The first inequality follows from the positive integer \(b_c-a_c\ge1\).
Thus the reduced rational core itself, not just a removable matrix scale,
has unbounded height as the four-point arc shrinks. Both output-conformal
normalization and actual source reanchoring preserve \(r=N\beta/N\alpha\).

Anchor the first phase at one. If its angular diameter is \(\Delta\le2\),
then \(t_0=0\), \(t_3=\tan(\Delta/2)\), and
\(\varepsilon\le\Delta\). For an endpoint cluster with squared radius
\(N\) and \(\Delta\le C N^{-1/4}\), \(C\le2\), (3) gives
\[
 b_c\ge C^{-2}\sqrt N.
\]
The target least squared radius nevertheless equals \(N\), because the
source and target finite sets coincide.

## 3. Actual integral endpoint examples

The primitive common-unit Pell family in
[four_point_bonus_counterexample.md](four_point_bonus_counterexample.md)
has four actual Gaussian lattice points with shrinking normalized arc
constant. It supplies infinitely many such source sets with unbounded
radius and eventually \(C<2\). Applying (1) to each set therefore gives
actual source-and-target endpoint examples of arbitrarily high coefficient
core and exactly unchanged radius. This assertion concerns four points,
not an unbounded cluster.

An independent exact rational calculation applied the construction to the
Pell indices \(1,11,21\), using the existing construction in
[check_gaussian_near_real_divisor_reduction.py](check_gaussian_near_real_divisor_reduction.py).
All four interchange identities, positive determinant, nonconformality, and
bounds (2)--(3) passed. The squared source radii had bit lengths
\(29,279,529\); the reduced core denominators had bit lengths
\(19,185,352\), respectively. Radius preservation is exact set equality,
not a numerical estimate. These cases are retained in
[the persistent exact checker](check_mobius_four_point_selfmap.py).

## 4. What the graph equation does and does not add

After choosing corresponding source and target anchors, a fractional-linear
map of affine complex circle coordinates has an exact relation
\[
 AXY+BX+CY=0,
\]
where \(X,Y\) are displacements from the anchors. Equivalently,
\[
 (AX+C)(AY+B)=BC.
\]
These equations are useful integrality constraints when their coefficients
are controlled. They also hold identically for the self-maps above, with
both sets contained in actual endpoint arcs and with an unbounded invariant
core. Thus neither this product identity nor exact four-point rational
interpolation can by itself force radius growth for all high-height cores.
The obstruction involves the pole between sampled points: although the
finite target set is short, the image of the whole source interval wraps
through the complementary arc. An interval-wide small-derivative premise
would exclude this phenomenon, but finite endpoint hypotheses do not imply
that premise.

A result for five or more points would need to exploit additional
correspondences or exclude such finite-set projective symmetries by a new
arithmetic argument. No such many-point conclusion is proved in this note.
