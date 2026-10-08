# A quadratic boundary barrier for eight-point conic Gale data

This note rules out one natural quadratic certificate for the rank-three
eight-point relation lattice.  Work over a field of characteristic zero.
Let `C` be the locus in `Gr(3,8)` whose eight row points lie on a conic,
including arbitrary nonzero row scalings.  For each unordered balanced
partition `S|S^c`, let `B_S` be the reducible-conic boundary on which the
four rows in `S` lie on one line and the four rows in `S^c` lie on another.
Then, in fact,

```text
I(union_(|S|=4) B_S)_2 = I(C)_2 = I(Gr(3,8))_2.     (1)
```

Thus a homogeneous quadratic in the triple Plucker coordinates which
vanishes on all 35 balanced boundaries already vanishes on every smooth
conic configuration.  In particular there is no nonzero quadratic class,
modulo the Plucker and conic relations, having the proposed boundary
vanishing property.

Section 4 then imposes both arithmetic constraints omitted from this first
statement: the Gale rows sum to zero and their scalings have the form
`q_i^2/D_i`.  The same quadratic obstruction survives.

## 1. Row-weight decomposition

Write a smooth conic row as

```text
g_i=a_i(s_i^2,s_i t_i,t_i^2),
[ij]=s_i t_j-t_i s_j.
```

For an increasing triple `I={i,j,k}`, its Plucker coordinate is

```text
p_I=a_i a_j a_k [ij][ik][jk].                       (2)
```

The diagonal torus which independently changes the eight `a_i` decomposes
the quadratic coordinate ring into row-weight spaces.  A monomial
`p_I p_J` has weight

```text
m=(1_I(1)+1_J(1),...,1_I(8)+1_J(8)).
```

If a quadratic vanishes on every `B_S`, varying the nonzero row scalings on
that boundary shows that each of its row-weight components vanishes there
separately.  It is therefore enough to compare boundary and smooth-conic
evaluation spans in each weight space.

Put `r=|I intersect J|`.  Up to relabeling, the weight spaces have the
following sizes:

| `r` | number of monomials in the weight space | number of weights |
|---:|---:|---:|
| `0` | `10` | `28` |
| `1` | `3` | `280` |
| `2` or `3` | `1` | `476` |

There are `1596` monomials in all, as expected.

## 2. The smooth-conic ranks are at most `1,2,5`

The one-dimensional cases require no argument.  In the `r=1` case, take
the repeated label to be `1` and the other labels to be `2,3,4,5`.  Order
the three monomials as

```text
p_123 p_145,  p_124 p_135,  p_125 p_134.             (3)
```

After (2), all three share the factor
`a_1^2 a_2 a_3 a_4 a_5 product_(j=2)^5 [1j]`.
The remaining factors are the three matchings on labels `2,3,4,5`, so the
four-point Plucker relation gives

```text
m_0-m_1+m_2=0.                                      (4)
```

Hence the restriction rank is at most two.

For `r=0`, use labels `1,...,6` and order the ten complementary-triple
monomials as

```text
m_0=p_123 p_456,   m_1=p_124 p_356,
m_2=p_125 p_346,   m_3=p_126 p_345,
m_4=p_134 p_256,   m_5=p_135 p_246,
m_6=p_136 p_245,   m_7=p_145 p_236,
m_8=p_146 p_235,   m_9=p_156 p_234.                 (5)
```

Repeated use of `[ij][kl]-[ik][jl]+[il][jk]=0`, or direct homogeneous
expansion of (2), gives the five identities

```text
-m_0+m_1-m_2+m_3=0,
 m_0+m_4-m_5+m_6=0,
 m_0-m_1+m_2+m_4-m_5+m_7=0,
 m_0+m_2-m_5+m_8=0,
-m_0+m_1-m_4+m_9=0.                                 (6)
```

They are independent because `m_3,m_6,m_7,m_8,m_9` occur as distinct
pivots.  The smooth-conic restriction rank is therefore at most five.

## 3. Exact boundary witnesses attain those ranks

Parameterize a dense open set of `B_S` by taking row scaling one and

```text
g_i=(1,t_i,0)  for i in S,
g_i=(1,0,t_i)  for i outside S.                     (7)
```

Both groups lie on the two components of the reducible conic `YZ=0`.
All parameters below are nonzero and are pairwise distinct on each line.

For the `r=1` ordering (3), use `t=(1,2,3,4,5,6,7,8)`.  The cuts
`S=1236` and `S=1246` give evaluation rows, up to nonzero scalar,

```text
(0,1,1),
(1,0,-1).
```

The minor in the first two columns has determinant `-1`.  Thus the
boundary evaluation rank is two, equal to the smooth upper bound.

For the `r=0` ordering (5), the first four rows below use
`t=(1,2,3,4,5,6,7,8)` and the indicated cut.  The last uses
`t=(5,2,3,4,5,6,7,8)`.

```text
S=1234: (0, 0, 1, 1, 0, 4,  4,  3,  3,  0)
S=1235: (0, 1, 0,-1, 3, 0, -3, -2,  0,  2)
S=1236: (0,-3,-3, 0,-8,-8,  0,  0, -5, -5)
S=1245: (1, 0, 0, 1,-9,-8,  0,  0, -9, -8)
S=1234: (0, 0,-3,-3, 0,-4, -4, -1, -1,  0)
```

The minor in columns `m_0,m_1,m_2,m_4,m_5` has determinant `8`.
Consequently the boundary evaluation rank is at least five and hence is
exactly the smooth-conic rank.  Relabeling gives the same conclusion for
all 28 weights of this type.  In every one-dimensional weight, a balanced
cut can split both triples, and a generic choice in (7) makes the sole
monomial nonzero.

It follows that, weight by weight, the span of evaluations on the union of
the 35 balanced boundaries equals the span of evaluations on `C`.  Taking
annihilators proves (1).  The rank bounds came from the symbolic identities
(4) and (6); the exact samples are used only for the matching lower bounds,
so no inference of a polynomial identity from finite sampling is involved.

The accompanying checker
[`check_eight_point_gale_quadratic_boundary.py`](check_eight_point_gale_quadratic_boundary.py)
expands every bracket polynomial over the integers, verifies (4) and (6),
checks the displayed boundary rows and determinants, and enumerates all
`784` row weights.

## 4. The actual `q_i^2/D_i` boundary still spans every quadratic

The circle moment kernel lies in the hyperplane

```text
H={x_1+...+x_8=0},
```

because the constant row belongs to the moment matrix.  Equivalently, the
eight rows of any Gale basis satisfy `sum_i g_i=0`.  Impose also the actual
row-scaling formula.  Let `A_S` be the clean-cut boundary constructed
below.  Then

```text
I(union_(|S|=4) A_S)_2 = I(Gr(3,H))_2.               (8)
```

Choose nonzero, pairwise distinct parameters `x_i` for `i in S` and `y_j`
for `j outside S`.  Put

```text
r_i=x_i                    (i in S),
r_j=epsilon/y_j            (j outside S),
q_epsilon(r)=A r^2+B r+epsilon C.                    (9)
```

Use `v_i=(r_i,1)`, `D_i=product_(k!=i) det(v_i,v_k)`, and the actual Gale
rows

```text
(q_epsilon(r_i)^2/D_i)(r_i^2,r_i,1).
```

For `epsilon!=0`, apply the common invertible coordinate change which sends
`(r^2,r,1)` to `(r,r^2,epsilon)`, and absorb `r_i` into the row scaling.
As `epsilon` tends to zero, the rows become

```text
b_i(1,x_i,0)     (i in S),
b_j(1,0,y_j)     (j outside S),                     (10)
```

where, up to one common sign,

```text
b_i=(A x_i+B)^2 /
    [x_i product_(k in S-{i})(x_k-x_i)],

b_j=(C+B/y_j)^2 /
    [y_j (product_(i in S)x_i)
     product_(k outside S, k!=j)(1/y_k-1/y_j)].      (11)
```

This is a formal algebraic degeneration, and below it is realized
p-adically. It is not a degeneration through positive real forms:
the limiting determinant is `-B^2/4`.

The Lagrange identities before specialization show that these eight rows
sum to zero.  Locally at `epsilon`, labels outside `S` have
`v(q_j)=v(epsilon)` and their mutual determinants have the same valuation;
all other `q_i` and determinants are units for generic parameters.  With
`B` a unit, the Gram determinant is also a unit.  Thus this is precisely a
clean balanced cut: `lambda_ij=v(epsilon)` on crossing edges and zero on
noncrossing edges.

Use the basis `e_1-e_8,...,e_7-e_8` of `H`, so the relevant Plucker
coordinates are the 35 minors on the first seven rows. The same row-weight
count now has 245 one-monomial weights, 105 three-monomial weights and
seven ten-monomial weights. The ordinary Grassmann--Plucker identities
(4) and (6), valid for arbitrary rank-three matrices, bound their ranks
by one, two and five. Therefore the quadratic evaluation rank is at most
`245+105*2+7*5=490`. Equivalently, the usual dimension formula gives

```text
dim S_(2,2,2)(k^7)=490,                             (12)
```

the dimension of the degree-two coordinate ring of `Gr(3,H)`.  The
deterministic checker
[`check_eight_point_actual_gale_quadratic_boundary.py`](check_eight_point_actual_gale_quadratic_boundary.py)
evaluates (10)-(11) over `F_1009`.  Its first 490 boundary points are
linearly independent among the 630 quadratic monomials in the 35
coordinates.  Every sample lifts to a rational point by using the same
integer parameters in (11), since every displayed denominator is nonzero
modulo 1009.  The determinant-square and positivity requirements can also
be retained in the lift.  Indeed `1009=1 mod 4`, so for prescribed nonzero
residues `A,B,C`, choose integer lifts with `A>0` and choose a positive
integer `S` solving the congruence

```text
S^2 = -B^2+4*1009*A*C mod 1009^2.
```

The nonzero square root of `-B^2` modulo 1009 lifts to this congruence
because its derivative is a unit. Define, exactly,
`C_tilde=(S^2+B^2)/(4*1009*A)` and `R=S/2`.
Then `C_tilde` is 1009-integral with residue `C`, and
`q(r)=A r^2+B r+1009 C_tilde` has determinant `R^2>0`.
The form is positive definite since `A>0`. This constructs rational
positive square-determinant Gram data with the specified boundary
reduction; no rational root of a fixed quadratic equation is being
inferred from Hensel lifting. Thus the modular minor also lifts to
actual circle Gram data. Independently, applying the same integer
parameters directly to the rational boundary formula (11) lifts the
nonzero modular determinant to a nonzero rational boundary determinant.
This is a certified lower bound of 490, while the preceding elementary
upper bound is 490, and proves (8).

The positive arithmetic lift works at every depth `epsilon=1009^e`.
For any desired further accuracy `k`, solve
`S^2=-B^2+4*1009^e*A*C mod 1009^(e+k)` and define
`C_tilde=(S^2+B^2)/(4*1009^e*A)`. Then `C_tilde=C mod 1009^k`.
With `k` tending to infinity along with `e`, the transformed actual
circle Gale rows converge p-adically to any of the displayed boundary
rows. Their mixed Plucker coordinates include a unit, so passing to
primitive integer Plucker coordinates changes their p-adic valuations
only by a common unit.

In particular a fixed homogeneous integral quadratic `F` for which
`N divides F(P)` on every actual primitive circle configuration would
vanish on these boundary rows: at depth `e` its valuation is at least
`e`, and passage to the p-adic limit forces zero. The rank-490 result
then makes `F` a Grassmann--Plucker relation on `Gr(3,H)`. The same
conclusion holds for a fixed loss `N divides c F(P)` with nonzero
integer `c`, since its finite 1009-adic valuation does not affect
the limit.

Therefore, even on the actual clean balanced-cut boundary, every
homogeneous quadratic which vanishes on all 35 cuts is already zero modulo
the Plucker relations and the linear equations expressing `K subset H`.
It cannot be nonzero on a smooth circle Gale configuration.  Unbalanced-cut
valuations and content need no further audit for this quadratic route,
because the required nonzero quadratic does not exist. This rules out
this homogeneous quadratic boundary-vanishing certificate. It does not
disprove a radius-versus-covolume inequality obtained by another argument.
