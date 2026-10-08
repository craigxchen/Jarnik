# Exact $M_{0,5}$ coordinates for the five-row CM cell

This note records the rational moduli coordinates used for the centers

\[
A=O,\qquad B=P,\qquad C=3P,\qquad D=7P
\]

on

\[
E:\ y^2=x^3-2x,\qquad P=(2,2).
\]

The two pencils of lines through `A` and `B` give the two coordinates
of the chart on \(\overline M_{0,5}\).  Normalize the first pencil by

\[
AC\longmapsto0,\qquad AD\longmapsto1,\qquad AB\longmapsto\infty,
\]

and the second by

\[
BC\longmapsto0,\qquad BD\longmapsto1,\qquad BA\longmapsto\infty.
\]

Write `x=x(Q)`, and let `s=s(P,Q)` be the slope of the line through
`P` and `Q`.  The exact constants are

\[
\begin{aligned}
x_P&=2, & x_C&=338,\\
x_D&=\frac{9374597220578}{5627138321281},\\
s_C&=\frac{1553}{84}, &
s_D&=\frac{927358740467713}{358778440454952}.
\end{aligned}
\]

Thus, away from the indicated poles,

\[
 a(Q)=\frac{(x-x_C)(x_D-x_P)}{(x-x_P)(x_D-x_C)},
 \qquad
 b(Q)=\frac{s(P,Q)-s_C}{s_D-s_C},
 \qquad
 s(P,Q)=\frac{y(Q)-2}{x(Q)-2}.
\]

The first formula is a Möbius transform of the `x`-coordinate, since
the pencil through `O` consists of vertical lines.  The second is a
Möbius transform of the chord slope through `P`.  At the basepoints,
use the tangent or pole limits:

\[
\begin{array}{c|cc}
Q&a(Q)&b(Q)\\ \hline
O&\displaystyle\frac{5594283994}{5632732605275}&\infty\\[6pt]
P&\infty&\displaystyle\frac{5736183875369054}{5705771236038721}\\[6pt]
-P&\infty&\infty.
\end{array}
\]

In particular `C` maps to `(a,b)=(0,0)`, and `D` maps to `(1,1)`.  The
line `AB` is contracted by the two pencils and its third intersection
with `E` is `-P`, explaining the simultaneous `(∞,∞)` value there.

## The ten boundary indices

A line through `rP` and `sP` has third intersection

\[
-(r+s)P,
\]

so the six lines through pairs of the four centers contribute

\[
-1,-3,-7,-4,-8,-10.
\]

Together with the four centers `0,1,3,7`, this gives the ten indices

\[
\mathcal K=\{0,1,3,7,-1,-3,-7,-4,-8,-10\}.
\]

Use the ordered directions

\[
(x_1,x_2,x_3,x_4,x_5)=(\infty,0,1,a,b).
\]

The exact evaluations and the resulting stable boundary are:

\[
\begin{array}{c|c|c|c}
k & (a(kP),b(kP)) & \text{colliding directions} & \text{boundary}\ \hline
0 & (\text{finite},\infty) & \{\infty,b\} & D_{15}\\
1 & (\infty,\text{finite}) & \{\infty,a\} & D_{14}\\
3 & (0,0) & \{0,a,b\}\text{ on the blown-up center} & D_{13}\\
7 & (1,1) & \{1,a,b\}\text{ on the blown-up center} & D_{12}\\
-1 & (\infty,\infty) & \{\infty,a,b\}\text{ on the blown-up center} & D_{23}\\
-3 & (0,\text{finite}) & \{0,a\} & D_{24}\\
-7 & (1,\text{finite}) & \{1,a\} & D_{34}\\
-4 & (\text{finite},0) & \{0,b\} & D_{25}\\
-8 & (\text{finite},1) & \{1,b\} & D_{35}\\
-10 & (a=b) & \{a,b\} & D_{45}
\end{array}
\]

At each of the three blown-up centers, the displayed triple is the
colliding cluster; the boundary label is the complementary two-element
subset.  For example, at `(0,0)` the cluster is `{2,4,5}`, so the
stable boundary is `D_{13}`.

The exact checker
[`check_five_row_m0_5_coordinates.py`](./check_five_row_m0_5_coordinates.py)
computes the group law with `Fraction`, checks the normalizations, checks
that the ten multiples are distinct, evaluates the ten boundary cases,
and tests sample multiples (N=2,4,5,6,8,9,11) away from the boundary.
It is a finite coordinate verification; it does not establish a uniform
count bound or any height estimate.
