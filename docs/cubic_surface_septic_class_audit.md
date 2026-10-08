# Cubic-surface class audit for the degree-seven septic

Let `S=Bl_(A_1,...,A_6)(P^2)` be a smooth cubic surface, with `L` the
pullback of a line and `E_i` the exceptional classes.  The anticanonical
hyperplane class is

\[
H=3L-E_1-\cdots-E_6.
\]

For a divisor `D=aL-sum b_i E_i`, a smooth genus-one curve of degree
seven in the cubic-surface embedding must satisfy

\[
H\cdot D=3a-sum b_i=7,
\qquad
D^2=a^2-sum b_i^2=7.
\]

The second equation is adjunction: `2g-2=D.(D+K_S)=D^2-H.D=0`.
Cauchy gives

\[
(3a-7)^2\le 6(a^2-7),
\]
so `3 <= a <= 11`.  Exact enumeration of the resulting finite lattice
problem gives several integral coefficient patterns; therefore the class
and genus equations alone do not obstruct a cubic-surface septic.

There is also a concrete smooth family.  Choose a smooth plane cubic
`F=0` through `A_5,A_6`, avoiding `A_1,...,A_4`, with transverse
incidence at the two blown-up points.  Its strict transform has

\[
D=3L-E_5-E_6,
\quad H\cdot D=7,
\quad D^2=7,
\quad g(D)=1.
\]

The strict transform is smooth, and the anticanonical map embeds `S` as
a smooth cubic surface in `P^3`; hence it gives an embedded
genus-one degree-seven curve on a smooth cubic surface.

For a fully rational instance, take `E: y^2 z=x^3-x z^2`,
`A_5=(0:0:1)`, `A_6=(1:0:1)`, and
`A_1=(-3:-3:1), A_2=(-3:-2:1), A_3=(-2:-3:1), A_4=(-2:1:1)`.
These six points have no three collinear and do not lie on a conic, so
their blowup is a smooth cubic surface.  The strict transform of E is the
class above.  This is an explicit septic on a smooth cubic, but it is not
an example of the balanced six-point trisecant system.

This construction does not satisfy the ten balanced trisecant conditions
automatically.  For five selected points (p_i) on the curve, each line
(overline{p_ip_j}) has a residual third intersection with the cubic
surface, and requiring that residual point to lie on the curve remains an
additional incidence condition.  No existence claim for all ten residual
contacts follows from the divisor calculation.
