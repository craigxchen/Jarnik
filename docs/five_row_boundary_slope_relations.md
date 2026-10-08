# Exact five-row boundary slope table in one anticanonical section

This note fixes one common `SL_2` chart and computes all ten boundary
restriction pairs for one section `Q`.  It is an algebraic supplement to
[five_row_boundary_node_height.md](five_row_boundary_node_height.md); it
does not introduce a second section or a numerical lattice configuration.

Write
`Delta_ij=Y_i X_j-X_i Y_j`.  Thus on an outside pair `(X_i,Y_i)=(1,z_i)`
it is `z_i-z_j`, and an inside row at `(0,1)` against an outside row at
`(1,z_j)` contributes `+1`.  Use the following six independent
degree-two-each graph coordinates (indices here are one-based):

```text
G0 = Delta12 Delta23 Delta31 Delta45^2
G1 = Delta12 Delta24 Delta41 Delta35^2
G2 = Delta12 Delta25 Delta51 Delta34^2
G3 = Delta13 Delta34 Delta41 Delta25^2
G4 = Delta13 Delta35 Delta51 Delta24^2
G5 = Delta14 Delta45 Delta51 Delta23^2.
```

Define one section by `Q=sum(q_r G_r)`, with
`q=(q0,q1,q2,q3,q4,q5)`.  For an inside pair `S`, order its outside
labels increasingly as `(a,b,c)`, put those two rows at `(0,1)`, and write
the restriction in the common affine chart as

```text
Q|_(D_S) = alpha_S (z_a-z_c) + beta_S (z_b-z_c).
```

The exact linear forms are:

| `S` | outside `(a,b,c)` | `alpha_S` | `beta_S` |
|---|---|---|---|
|12|345|`-q3-q4`|`q3-q5`|
|13|245|`-q1-q2`|`q1-q5`|
|14|235|`-q0-q2`|`q0-q4`|
|15|234|`-q0-q1`|`q0-q3`|
|23|145|`q1+q2+q3+q4`|`-q1-q3`|
|24|135|`q0+q2-q3+q5`|`-q0+q3`|
|25|134|`q0+q1-q4-q5`|`-q0+q4`|
|34|125|`-q0-q1+q4+q5`|`q0+q1`|
|35|124|`-q0-q2+q3-q5`|`q0+q2`|
|45|123|`-q1-q2-q3-q4`|`q1+q2`|

The table is obtained by exact specialization of the six graph monomials;
the checker verifies both the restriction identity and that these six
coordinates have rank six.

## Shared-node identities

If disjoint pairs `S,T` meet at a node and `c` is the remaining label, the
two restrictions use the same `Q` and the same chart.  The coefficient of
`z_c` on the two sides is opposite.  Translating that statement through the
table gives the fifteen exact identities

```text
(alpha12+beta12)+(alpha34+beta34) = 0
(beta12)-(alpha35+beta35)          = 0
(alpha12)-(alpha45+beta45)         = 0
(alpha13+beta13)+(alpha24+beta24) = 0
 beta13-(alpha25+beta25)          = 0
 alpha13+beta45                   = 0
(alpha14+beta14)+(alpha23+beta23) = 0
 beta14+beta25                    = 0
 alpha14+beta35                   = 0
 beta23-(alpha15+beta15)            = 0
 beta15+beta24                    = 0
 alpha15+beta34                   = 0
 alpha23+alpha45                   = 0
 alpha24+alpha35                   = 0
 alpha25+alpha34                   = 0.
```

For example, at `D_12 intersect D_35`, the remaining label is `4`.
It is the second outside coordinate on `D_12`, so its coefficient is
`beta12`; it is the third outside coordinate on `D_35`, so its coefficient
is `-(alpha35+beta35)`.  In the fixed graph trivializations these opposite
node coefficients sum to zero, giving the second identity.

These are the exact degree-one compatibilities available from one `Q`.
For each of the fifteen disjoint pairs `(S,T)`, the four linear forms
`(alpha_S,beta_S,alpha_T,beta_T)` have exact rank three.  The displayed
node identity is therefore their unique linear relation up to scale.  This
also proves the precise slope-ratio statement.  Let
`u_S(alpha_S,beta_S)` and `u_T(alpha_T,beta_T)` be the two node-coordinate
forms occurring in that identity.  The rank-three image is the hyperplane

```text
u_S + u_T = 0
```

in four-space.  Given any two generic projective slope pairs
`[alpha_S:beta_S]` and `[alpha_T:beta_T]` for which both node forms are
nonzero, choose nonzero scale factors `s,t` with
`s u_S+t u_T=0`.  The rank-three map then lifts
`[s alpha_S:s beta_S:t alpha_T:t beta_T]` to a section `Q`.  Hence the
map from the nonnodal open in the `Q`-space to
`P^1 times P^1` of the two slope ratios is dominant: no algebraic equation
in those two ratios is forced by their common node.  Any arithmetic coupling
of the two boundary congruences must retain their scales or use relations
involving additional boundary lines.

For a non-nodal section, every node coefficient in these identities is
nonzero.  The separate condition that the boundary zero is not a node is
`alpha_S+beta_S != 0`; the other two boundary points in the degree-one
restriction are `alpha_S != 0` and `beta_S != 0`, as in (8) of the boundary
height note.  These nonvanishings concern this one section `Q` throughout.

Run the exact standard-library verification with:

```text
python3 -B docs/check_five_row_boundary_slope_relations.py
```
