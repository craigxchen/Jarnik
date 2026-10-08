# The five-row boundary table as a Petersen circulation

The ten pair labels `S subset {1,...,5}`, `|S|=2`, are the vertices of the
Petersen graph.  Two vertices `S,T` are adjacent when they are disjoint; the
remaining label is `c`.  If

```text
F_S = alpha_S(z_a-z_c') + beta_S(z_b-z_c')
```

uses the increasing outside ordering `(a,b,c')`, define the edge flow on
`S--T` to be the coefficient of `z_c` in `F_S`.  The coefficient of `z_c`
in `F_T` is its negative by the common-chart node identity, so the flow is
well-defined after orienting each edge from the smaller to the larger vertex
index.

At a vertex `S`, its three incident node coefficients are
`alpha_S`, `beta_S`, and `-alpha_S-beta_S`.  Their sum is zero.  Therefore
the 15 edge flows form an integral circulation on the Petersen graph.

The exact checker establishes the following global facts:

* the oriented Petersen incidence matrix has rank `9`;
* every triangle-times-squared-edge graph gives a signed six-cycle flow, and
  every pentagon graph gives a signed five-cycle flow;
* the map from the six displayed `q`-coordinates to the 15 edge flows has
  rank `6`, and its gcd of nonzero `6 by 6` minors is `2`;
* using all 22 graph coordinates, the flow image contains a `6 by 6` minor
  equal to `-1`, so the complete graph flow lattice has index `1`;
* the 15 gluing equations among the 20 entries `(alpha_S,beta_S)` have rank
  `14`; and
* with all ten restrictions written in the same variables `z_1,...,z_5`,
  the exact identity `sum_S F_S(z)=0` holds coefficientwise.

Since the connected Petersen cycle lattice is the saturated integer kernel of
the incidence matrix, the gcd-of-minors computation says that the lattice of
flows from the six displayed graph coordinates has index exactly `2` in the
full integral cycle lattice.  This is a sublattice statement, not a claim
about every integral anticanonical section.  After all 22 graph coordinates
are allowed, the exhibited determinant `-1` gives the saturated cycle
lattice.  The separate unimodular coefficient-basis check for the integral
invariant space is recorded in
[five_row_gradient_discriminant_core_count.md](five_row_gradient_discriminant_core_count.md).
Modulo `2`, every flow from the first six coordinates has even total edge
support, since each is a six-cycle.  The twelve five-cycle generators among
the remaining graph coordinates have odd parity and remove this index-two
obstruction.

For the fixed section `q=(1,2,3,4,5,6)`, the ten common-variable linear
forms `F_S(z)` have rank `4`, the maximum allowed by translation invariance
(`sum_i` of each coefficient vector is zero).  This is a witness about one
section; it is not a genericity or arithmetic conclusion.

The [cycle equations](five_row_petersen_slope_cycles.md) eliminate the
vertex scales and give five necessary and sufficient slope compatibilities
on the nonnodal open set, including an explicit five-cycle relation.

The circulation and gluing identities are additive algebraic constraints.
They do not identify the distinct Gaussian moduli or their congruence
moduli, so no multiplication of boundary congruences or uniform arc bound is
claimed here.

Run the exact standard-library checker with:

```text
python3 -B docs/check_five_row_petersen_circulation.py
```
