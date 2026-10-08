# Balanced six-point elliptic curves as trisecant septics

This note verifies the Kapranov numerical class and gives an exact
linear-series formulation of the remaining existence question. It does
not establish existence or nonexistence of a balanced genus-one curve
in `bar(M)_(0,6)`.

## 1. The class and the meaning of the ten contacts

Use the Kapranov model obtained by blowing up five points
`P_1,...,P_5` in linear general position in `P^3`, then the strict
transforms of their ten connecting lines. Denote the pullback
hyperplane by `H`, the point exceptional divisors by `E_i`, and the
line exceptional divisors by `E_ij`. This model and the plane class
below are recorded in Castravet, *The Cox ring of Mbar_0,6*,
[Notation 1.2 and the paragraph preceding it](https://arxiv.org/abs/0705.0070).
The twenty-five boundary divisors
are these fifteen exceptional divisors and the ten strict transforms
of the planes through triples of the points. The latter have classes

```text
B_ijk=H-E_i-E_j-E_k-E_ij-E_ik-E_jk.                    (1)
```

Let `f:E -> bar(M)_(0,6)` be birational onto its image, with `E` a
smooth genus-one curve and image meeting the open moduli space.
If every boundary pullback has degree one, then

```text
deg f^*E_i=deg f^*E_ij=1,
deg f^*H=7.                                          (2)
```

Indeed, any equation (1) gives `1=deg f^*H-3-3`.
The composite `phi:E->P^3` is therefore birational onto an integral
curve of degree seven. It is not presumed to be an embedding.

The point ideal of `P_i` pulls back to a divisor of degree one, so
there is exactly one point `p_i` of `E` above `P_i`, and the minimum
vanishing order of the three local coordinates is one. The branch is
immersed there. The line ideal of `P_iP_j` has total transform

```text
I_(P_iP_j) O_X=O_X(-E_i-E_j-E_ij).                    (3)
```

To see the endpoint contribution locally, write the line as
`y=z=0` and an immersed branch through the endpoint as
`x=t+..., y=t^a u(t), z=t^b v(t)`, after a suitable coordinate
choice when the tangent is the line. Blowing up the point changes
the ideal order from `min(a,b)` to `min(a,b)-1` along its strict
transform. Away from the two endpoints no point correction occurs.
The subsequent blow-up of that strict line accounts for the remaining
order. The ten strict lines are disjoint after blowing up the points.

Consequently the complete inverse-image divisor of each original
connecting line is

```text
phi^*I_(P_iP_j)=O_E(-p_i-p_j-q_ij)                    (4)
```

for one point `q_ij`, counted with multiplicity. Thus each connecting
line is a trisecant in the scheme-theoretic sense, of exact contact
length three. The point `q_ij` can equal an endpoint, in which case
the extra contact is tangency at that endpoint. It is not legitimate
to assume ten new distinct contact points without an additional
transversality hypothesis.

If every `q_ij` avoids the five centers, then those ten extra points
are distinct: lines with a common vertex only intersect at that
vertex, and connecting lines on four different vertices are skew
because the centers are in linear general position.

Each plane through three centers has degree-seven pullback containing

```text
p_i+p_j+p_k+q_ij+q_ik+q_jk                            (5)
```

as a divisor, with multiplicities interpreted through (3). Its final
strict boundary transform contributes the remaining degree one.
This statement does not require the six displayed points to be
distinct.

Conversely, suppose a basepoint-free map `phi:E->P^3`, birational onto
its image, has degree seven; its image contains the five prescribed
centers, with point-ideal pullbacks of degree one; and all ten
line-ideal pullbacks have degree three. Suppose its image is contained
in none of the ten boundary planes. The unique lift through the point
and line blow-ups then has degree one against every `E_i,E_ij`, and
(1) gives degree one against every `B_ijk`. These conditions are
therefore an exact existence criterion, not just necessary degree
counts.

## 2. An elementary restriction: no containing quadric

Such a septic lies on no quadric surface, including a reducible or
nonreduced quadric. If a quadric contained the curve, its restriction
to each connecting line would vanish on the length-three divisor in
(4); hence that degree-two restriction would vanish identically.
It would contain all ten connecting lines. In fact the six lines
between four linearly independent centers already force its defining
quadratic polynomial to be zero: polarization vanishes on every pair
of these basis vectors, including the diagonal pairs.

Equivalently, for `L=phi^*O(1)` and its four-dimensional defining
subspace `V`, the multiplication map

```text
Sym^2(V) -> H^0(E,L^2)                                (6)
```

is injective. The dimensions are ten and fourteen. This is compatible
with a genus-one septic; it supplies no contradiction by itself.
The analogous argument does not exclude a cubic surface: three zeros
on a line are compatible with a nonzero cubic restriction.

## 3. An exact elliptic linear-series incidence problem

The existence question can be posed without a spatial embedding
assumption. Choose:

```text
E a smooth genus-one curve,
L in Pic^7(E),
V subset H^0(E,L),     dim V=4,
p_1,...,p_5 distinct points of E.
```

Require `V` to be basepoint free, with the five images a projective
frame in `P^3` (every four linearly independent). Birationality then
follows automatically: the mapping degree divides seven, and degree
seven would give a line as image, contradicting the independence of
the four sections. This does not imply that the map is an embedding.
For each pair define the two-dimensional pencil

```text
W_ij=V intersect H^0(E,L(-p_i-p_j)).                   (7)
```

The projective-frame condition ensures the dimension is exactly two.
The ten trisecant conditions are precisely that the fixed divisor of
each pencil `W_ij` have degree exactly three. Equivalently, there is
a point `q_ij` with

```text
dim(V intersect H^0(E,L(-p_i-p_j-q_ij)))=2,            (8)
```

and the pencil has no additional common zero beyond that divisor.
The equality is meant scheme-theoretically when `q_ij` equals an
endpoint. One must also require the point-ideal condition at each
`p_i` to have order one; for an embedded or everywhere immersed map
with distinct fibers over the centers this is automatic.

A degree-three divisor on an elliptic curve cuts out a subspace of
codimension three in `H^0(E,L)`, which has dimension seven. Thus (8)
is a concrete rank-at-most-two condition for a `4 by 3` evaluation
matrix, for each of the ten pairs. It can be written in an elliptic
function basis and checked by minors and their first derivatives for
repeated points. It is a finite incidence system whose solutions,
subject to the stated open conditions, are equivalent to the desired
balanced six-curves.

For fixed `E,L,p_i,q_ij`, one such rank condition has expected
codimension two in `Gr(4,7)`; allowing `q_ij` to move leaves an
expected single condition per pair. These ten conditions are highly
correlated. Their expected dimension is not a proof that they can
be solved simultaneously, nor that an intersection lies in the
required open set. In particular neither a naive Hilbert dimension
count nor the remaining singularity budgets of the forgetful images
settles this incidence problem.

## 4. A four-section interpolation form

There is a smaller exact system suited to symbolic or finite-field
experiments.  Put four centers at the coordinate points

```text
P_1=(1:0:0:0), ..., P_4=(0:0:0:1)
```

and the fifth at `P_5=(1:1:1:1)`.  Write a basis of `V` as
`s_1,...,s_4`.  For `{i,j,k,l}={1,2,3,4}`, the coordinate section
`s_i` must have divisor

```text
div(s_i)=p_j+p_k+p_l+q_jk+q_jl+q_kl+r_jkl,            (9)
```

where `r_jkl` is the degree-one contact with the strict transform of
the plane through `P_j,P_k,P_l`.

Conversely, start with `E`, a degree-seven line bundle `L`, distinct
points `p_1,...,p_5`, and six points `q_ij` for
`1<=i<j<=4`.  For each `i`, the first six terms on the right of (9)
have degree six.  On an elliptic curve the degree-one line bundle
obtained after subtracting them from `L` has a unique effective
divisor.  It therefore defines a unique section `s_i` up to scalar and
the residual point `r_jkl`.  Provided `s_i(p_5)` is nonzero, scale the
four sections so that

```text
s_1(p_5)=s_2(p_5)=s_3(p_5)=s_4(p_5).
```

The six connecting lines among `P_1,...,P_4` then already have the
required common-zero divisors

```text
gcd(div(s_k),div(s_l))=p_i+p_j+q_ij,                  (10)
```

where `{k,l}` is the complement of `{i,j}`.  Equality in (10), rather
than mere containment, is an open no-extra-common-zero condition.

Only the four lines from `P_5` remain.  For each `i`, take any two
independent differences among the three sections `s_j` with `j!=i`.
Their common divisor must be exactly

```text
p_5+p_i+q_i5.                                         (11)
```

They already vanish at `p_5` by the normalization and at `p_i` by
(9).  Thus (11) asks for exactly one further common zero, counted
scheme-theoretically, and no fourth one.  This gives four explicit
resultant or evaluation-minor equations, one for each `i`.

The ten plane divisors are then visible without further conditions.
Besides (9), for `{i,j,k,l}={1,2,3,4}` one has

```text
div(s_i-s_j)
 =p_5+p_k+p_l+q_5k+q_5l+q_kl+r_5kl.                  (12)
```

The last point is the residual zero of the degree-seven section.
Equations (9)--(12), together with basepoint freeness, simple point
ideals at the `p_i`, linear independence, and the no-extra-zero clauses,
are equivalent to the criterion in Section 3.  They show that the
unresolved question is a fourfold elliptic common-zero interpolation
problem.  The construction has substantial moving data, so the class
calculation and parameter counts alone give no nonexistence result.
