# Five-row invariant heights and the degree-five del Pezzo surface

The five-row invariant space has an exact boundary-height interpretation.
Its primitive projective height is `10w+o(w)`, while each of the ten
boundary divisor heights is `2w+o(w)`. A boundary divisor combines the
two complementary core cuts, and does not remember their separate
Gaussian orientations. The degree-one versus degree-five comparison
on an anticanonical elliptic section therefore matches the fair
profile; it does not itself yield a saving.

The particular irreducible five-cycle relation constructed earlier is
a special smooth elliptic section whose ten boundary intersections
are paired. It is excluded from an asymptotic full profile by its
two-monomial equation, before any general elliptic-height theorem is
used. The checker is
[check_five_row_del_pezzo_arithmetic.py](check_five_row_del_pezzo_arithmetic.py).

## 1. The geometric space and its line bundles

[Bauer--Catanese, Sections 1--2](https://link.springer.com/article/10.1007/s12215-020-00483-9)
identify the moduli space of stable ordered quintuples on the projective
line with the split del Pezzo surface of degree five. Its anticanonical
space has dimension six. Its ten boundary lines are indexed by pairs;
two distinct lines intersect precisely when those pairs are disjoint.
The paper also identifies the bracket pentagons with anticanonical
hyperplane sections. These give the geometric interpretation used
below; the arithmetic gcd calculations are explicit here.

Write this surface as `Y`, its anticanonical bundle as `L=-K_Y`, and
the boundary line indexed by a pair `e` as `D_e`. On the chart
`(infinity,0,1,a,b)`, a useful compactification is

```text
Y=Bl_((0,0),(1,1),(infinity,infinity))(P^1 x P^1),
L=(2,2)-E_0-E_1-E_infinity.                             (1)
```

Its sections are bidegree-at-most `(2,2)` polynomials vanishing at
`(0,0)` and `(1,1)`, with coefficient of `a^2 b^2` zero. Thus this
space has dimension `9-3=6`, as does the five-row degree-two-each
invariant space `V_5`. The degree is `L^2=8-3=5`.

For reference, the complete Hilbert function has the exact invariant
calculation

```text
dim V_(5,2d)=binomial(5d+3,3)-5binomial(3d+2,3)
              +10binomial(d+1,3)
            =5d(d+1)/2+1.                              (2)
```

It follows by multiplying the invariant weight-generating function
by `1-x` and using inclusion-exclusion; binomial coefficients with
upper entry below three are zero.

Each `D_e` satisfies `L.D_e=1` and `D_e^2=-1`. The line bundle
`L-D_e` is globally generated with four independent sections. For
the diagonal `D_45` in (1), its class is
`(1,1)-E_0-E_1-E_infinity`, so `L-D_45=(1,1)` is the pullback of
the globally generated system on `P^1 x P^1`. Permuting the five
labels proves the assertion for every boundary line.

## 2. Exact arithmetic of all twenty-two graph coordinates

Use the actual five-row central factorization, retaining its inherited
normalization, and put

```text
|Delta_ij|=b_ij product_(S containing i,j) n_S,
(1-eta)w<=log n_S<=(1+eta)w,
1<=b_ij<=exp(beta),       B=product_(10 edges)b_ij.      (3)
```

The empty cut may be omitted; it has zero exponent in every formula
below. The twenty-two degree-two-each bracket graph monomials are
the twelve Hamilton cycles and the ten triangle products times a
squared complementary edge. All are nonzero on distinct directions,
and they span `V_5`. Denote their numerical values by `v_G`, and put

```text
g=gcd_Graphs |v_G|,
D_0=product_S n_S^g_0(|S|),       g_0(a)=max(0,2a-5).
```

At every cut the minimum number of internal graph edges is `g_0`.
Each edge occurs at most twice in a graph. Choosing a graph attaining
the minimum separately at each prime therefore gives the exact
divisibilities

```text
D_0 | g,       g | D_0 B^2.                             (4)
```

The sum of the core exponents is

```text
sum_(a=0..5) binomial(5,a) g_0(a)=30.
```

Every graph has five brackets, each with eight incident core blocks,
so its raw logarithmic size is between `40(1-eta)w` and
`40(1+eta)w+5beta`. The primitive projective height `h_L` using these
twenty-two coordinates consequently satisfies

```text
10w-70eta w-20beta <= h_L <= 10w+70eta w+5beta.          (5)
```

Using a fixed six-coordinate basis changes this only by a bounded
constant. No cancellation of a sum is used in (4) or (5).

Fix an edge `e`. The graph monomials containing `e` span the
four-dimensional space of anticanonical sections vanishing on `D_e`.
Let `g_e` be the gcd of their numerical values. Their exact cut
minimum is

```text
g_e(S)=g_0(|S|)+1_(S=e or S=e^c).                       (6)
```

This minimum is attained for all thirty-two cuts, as the exact
checker verifies. In particular

```text
D_0 n_e n_(e^c) | g_e,
g_e | D_0 n_e n_(e^c) B^2,
g | g_e.                                               (7)
```

The positive integer `q_e=g_e/g` is the finite-place contact content
of the projective point with the boundary line, up to fixed choices
of an integral basis. Equations (4) and (7) give

```text
|log q_e-log n_e-log n_(e^c)| <= 2log B <=20beta.         (8)
```

This is the two-cut arithmetic dictionary. It retains prime powers
and explicitly charges extra gcds from corrections and primitive
residues. It is not a statement about one Gaussian orientation alone.

## 3. The archimedean term is also accounted for

The subspace in (6), after its common divisor `D_e` is removed,
is the basepoint-free system `L-D_e`. Let `h_(L-D_e)` be the
projective height of its graph-coordinate values. Define the
associated Weil height for `D_e` by their difference. It is exactly

```text
h_(D_e)=h_L-h_(L-D_e)
       =log q_e+log(max_all |v_G| / max_(G containing e)|v_G|).
                                                               (9)
```

The last term is nonnegative. Both maxima have a nonzero graph
monomial with the raw bounds used in (5), so it is at most
`80eta w+5beta`. Therefore

```text
2w-2eta w-20beta <= h_(D_e)
                    <=2w+82eta w+25beta.                (10)
```

Thus the boundary height `2w+o(w)` includes the archimedean place;
it is not merely the finite gcd in (8). There is also the geometric
identity `sum_e D_e ~ 2L`, consistent with the sum `20w+o(w)` of
the ten heights and `2h_L=20w+o(w)`.

For a smooth anticanonical section `E`, adjunction gives genus one,
`deg(L|E)=5`, and `deg(D_e|E)=1`. If a boundary line is not a
component, its restriction is one rational point counted with
multiplicity. Thus the expected ratio is exactly

```text
h_(D_e)/h_L -> 1/5,
```

which agrees with (5) and (10). For a fixed smooth section, standard
elliptic canonical-height comparison makes the difference from this
ratio a degree-zero height, generally of order `sqrt(h_L)` rather
than bounded. For varying section coefficients, a fixed-curve
constant cannot silently be made uniform; a separate bound in the
coefficient height is necessary. The arithmetic bounds (5), (8),
and (10) themselves are independent of any chosen section equation.
The exact comparison with those varying constants retained, and a
nonspecial fixed elliptic section with infinitely many rational points
realizing all ten aggregate proportions, are proved in
[five_row_elliptic_height_compatibility.md](five_row_elliptic_height_compatibility.md).

## 4. The earlier five-cycle relation is a special smooth pencil

Let

```text
C_1=Delta_12 Delta_23 Delta_34 Delta_45 Delta_51,
C_2=Delta_13 Delta_35 Delta_52 Delta_24 Delta_41,
G_t=t C_1+C_2.
```

On the chart in (1),

```text
G_t(a,b)=-t a^2+[t+(t+1)b-b^2]a-tb.                    (11)
```

For `t=12` this is the irreducible numerical relation from
[minimal_support_invariant_descent.md](minimal_support_invariant_descent.md).
It is genuinely smooth of genus one, not merely a generic assertion.
Its double-cover discriminant is

```text
b^4-26b^3+145b^2-264b+144
  =(b-1)(b-4)(b^2-21b+36),
disc=700710912 !=0.                                    (12)
```

The four branch points are distinct, so the normalization has genus
one. Since the integral anticanonical section has arithmetic genus
one, it has no singularity contributions and is smooth.

There is also an exact Tate form. With
`D=b-t(a-1)`, the map

```text
x=t b/D,       y=t^2/D
```

gives

```text
E_t: y^2+(1-t)xy-t y=x^3-tx^2,
Delta(E_t)=t^5(t^2-11t-1).                              (13)
```

The checker verifies the homogeneous identity, including all factors
of `t`. Its five boundary points are

```text
O, (0,0), (t,t^2), (t,0), (0,t),                        (14)
```

the multiples of the point `T=(0,0)` of order five. This can be
checked without a torsion classification: the tangent `y=0` at
`T` gives `2T=(t,t^2)`, while `y=tx` passes through `T` and is
tangent at `2T`, giving `5T=0`. The five displayed points are
distinct when `t!=0` and the discriminant is nonzero.

The geometric reason there are five points, rather than ten, is
that `C_1` and `C_2` are complementary boundary pentagons. Every
line in either pentagon meets exactly one line in the other,
giving five transverse boundary-node basepoints. Every smooth
member of their pencil passes through these basepoints. Since its
intersection with each boundary line has degree one, the two
boundary restrictions at each node coincide. This is a special
feature of the binomial pencil, not a property of a general
anticanonical section.

Explicitly, in homogeneous plane coordinates `(a:b:z)`, the paired
boundary points and their Tate images are

| Boundary pair | Plane point | Tate image |
|---|---|---|
| `D_12,D_35` | `(1:1:1)` | `2T=(t,t^2)` |
| `D_23,D_14` | `(1:0:0)` | `T=(0,0)` |
| `D_34,D_25` | `(1:0:1)` | `O` |
| `D_45,D_13` | `(0:0:1)` | `4T=(0,t)` |
| `D_15,D_24` | `(0:1:0)` | `3T=(t,0)` |

The homogeneous transformation is
`(x:y:z_new)=(t b:t^2 z:b-ta+tz)`; the checker verifies all these
projective images as polynomial identities.

The point `(a,b)=(2,3)` on `G_12` maps to `(-4,-16)` on `E_12`.
It is a point of order two, as the exact negation formula shows;
it should not be used to assert a positive-rank point sequence.
No claim about the full Mordell--Weil rank is made here.

## 5. The fixed binomial is already excluded by full-profile arithmetic

The two pentagon monomials have no common edge, and at every cut
their minimum internal-edge count is exactly `g_0`. Their numerical
gcd is therefore between `D_0` and `D_0 B`, where the upper
correction needs only the first power of `B`. Thus

```text
10w-70eta w-10beta <= h(C_2/C_1)
                       <=10w+70eta w+5beta.             (15)
```

If `G_t=0`, the ratio equals `-t`. Consequently the rational
coefficient height must satisfy

```text
h(t) >= 10w-70eta w-10beta.                             (16)
```

At an individual pair cut, the cycle containing that edge has one
extra core norm power; the same is true at the complementary
triple cut. Their equality therefore already requires those powers
in coefficients or correction/residue factors. Formula (16) is
the simultaneous projective-height version of that observation.
In particular fixed `t=12`, or any `h(t)=o(w)`, cannot produce the
full profile with negligible correction heights.

This exclusion is the earlier multiplicative projective-unit
obstruction expressed geometrically. It does not exclude a general
irreducible anticanonical equation with several graph monomials.
For such a curve the ten boundary points can be distinct, the
degree ratios match (10), and the moduli heights forget which of
the two complementary oriented Gaussian blocks supplied each
boundary contact. That lost orientation, together with the
actual primitive pair residues and the arithmetic chart, is the
remaining information needed for any transfer back to the
uniform lattice-arc problem.
