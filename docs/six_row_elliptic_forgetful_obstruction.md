# A forgetful obstruction for balanced elliptic six-curves

The five-row boundary-height ratio is compatible with a fixed smooth
anticanonical elliptic section.  At six rows, overlapping five-point
projections impose an additional constraint.  This note proves a
conditional geometric obstruction: an elliptic curve in
`bar(M)_(0,6)` whose twenty-five boundary divisors all have the same
positive degree can have at most two smooth forgetful images in
`bar(M)_(0,5)`.  At common degree one, at most three forgetful images
can even be anticanonical; this includes singular rational
anticanonical images.  If there are three, exactly one is smooth
elliptic and two are singular rational.

This does not assert that an arbitrary full independent profile lies on
a fixed elliptic curve, or even that it has one five-row relation.  It
classifies what would happen after such a fixed six-row elliptic support
has been justified.

## 1. Equal boundary heights become equal intersection degrees

Work over a characteristic-zero field and pass to its algebraic
closure when discussing deck transformations.  Let `C` be a smooth
genus-one curve and

```text
f:C -> bar(M)_(0,6)
```

a morphism birational onto an integral curve meeting the open moduli
space.  Suppose

```text
deg f^*D_I=d>0                                           (1)
```

for every one of the twenty-five stable boundary divisors, with the
same integer `d`.

This is the exact fixed-curve consequence of a balanced boundary-height
profile.  If `Q_N=NP+Q_0` runs through a non-torsion elliptic orbit,
then each fixed divisor has

```text
h_(D_I)(f(Q_N))=deg(f^*D_I) hat(h)(P)N^2+O(N).
```

Thus equality of all quadratic leading coefficients forces (1).  The
existence of `C,f` is an additional hypothesis, not a consequence of
the height asymptotic alone.

## 2. Every forgetful image has a forced divisor class

Fix a label `j` and let

```text
pi_j:bar(M)_(0,6) -> bar(M)_(0,5)
```

forget that label.  A boundary line downstairs is indexed by a
two-element subset `e` of the other five labels.  Its pullback is

```text
pi_j^*D_e=D_e+D_(e union {j}).                           (2)
```

Both terms are stable six-point boundary divisors.  Equation (1)
therefore gives

```text
deg (pi_j f)^*D_e=2d                                    (3)
```

for all ten boundary lines downstairs.

Let `Gamma_j` be the integral image curve, let `tilde(Gamma)_j` be its
normalization, and let `r_j` denote the degree of the induced map
`C->tilde(Gamma)_j`.  On the degree-five del Pezzo model

```text
Y=Bl_(p_1,...,p_4)(P^2),       L=-K_Y=3H-E_1-...-E_4,
```

write the class of `Gamma_j` as `aH-sum m_iE_i`.  If its common
intersection with the ten boundary lines is `delta_j`, intersection
with `E_i` gives `m_i=delta_j`, and intersection with every
`H-E_i-E_k` gives `a-2delta_j=delta_j`.  Hence, exactly,

```text
Gamma_j ~ delta_j L,       r_j delta_j=2d.              (4)
```

The arithmetic genus of this class is

```text
p_a(delta L)=1+5delta(delta-1)/2.                       (5)
```

If `Gamma_j` is smooth, Riemann--Hurwitz for the nonconstant map from
`C` shows that its genus is at most one.  Equations (4)--(5) then force

```text
delta_j=1,       Gamma_j in |L| smooth elliptic,
r_j=2d.                                               (6)
```

Thus every smooth forgetful image is an anticanonical elliptic curve,
independently of the value of the common boundary degree `d`.

Let `g_j` be the genus of the normalization of `Gamma_j`.
Riemann--Hurwitz gives `g_j in {0,1}`, but it does not force `g_j=1`.
The exact total singularity delta-invariant is

```text
p_a(Gamma_j)-g_j=1+5delta_j(delta_j-1)/2-g_j.           (7)
```

If `r_j=1`, the map from `C` to the normalization is birational and
therefore `g_j=1`.

For `d=1`, (4) gives the exact alternatives.  Either
`r_j=2,delta_j=1`, where the anticanonical image is smooth elliptic
(`g_j=1`, singularity budget zero) or singular rational (`g_j=0`,
budget one), or `r_j=1,delta_j=2`, where the image lies in `|2L|`, has
arithmetic genus six, and has genus-one normalization with total
singularity delta-invariant five.

## 3. Three smooth images are impossible

Normalize the six points on the open chart as

```text
(infinity,0,1,a,b,c).
```

Consider the three forgetful maps dropping respectively `a,b,c`.  Their
function fields are

```text
k(b,c),       k(a,c),       k(a,b).                     (8)
```

Suppose all three images are smooth.  By (6), each map is an unramified
isogeny of degree `2d` between genus-one curves.  Let its translation
kernel be `K_a,K_b,K_c`, respectively.  Each finite group has even
order and therefore contains a nonzero two-torsion point; choose

```text
t_a in K_a[2],       t_b in K_b[2],       t_c in K_c[2].
```

There is always a nonzero element in

```text
(K_b+K_c) intersect (K_a+K_c) intersect (K_a+K_b).      (9)
```

Indeed, if two of `t_a,t_b,t_c` agree, that common point lies in all
three sums.  If they are pairwise distinct, they are the three nonzero
points of `C[2]`, and every pair generates all of `C[2]`; then (9)
contains `C[2]`.

The function `a` is invariant under both `K_b` and `K_c`, hence under
their sum.  Similarly, `b` is invariant under `K_a+K_c`, and `c` under
`K_a+K_b`.  A nonzero translation from (9) therefore fixes all three
functions.  But `f` is birational onto its image, so

```text
k(C)=k(a,b,c),                                         (10)
```

and no nontrivial automorphism can fix that field.  This contradiction
proves that the three images cannot all be smooth.

The labels used as `a,b,c` were arbitrary.  If any three of the six
forgetful images were smooth, use the complementary three labels as the
normalizing anchors and repeat the argument.  Therefore

```text
at most two of the six forgetful images are smooth.     (11)
```

## 4. Degree one: at most four anticanonical images

Now assume `d=1`, and let `S` be the set of labels whose forgetful
images have `delta_j=1`.  These images may be smooth elliptic or
singular rational.  In either case (4) gives a degree-two map from `C`
to the normalization, hence a deck involution `tau_j` of `C`.

The involutions belonging to distinct labels in `S` are distinct.  If,
for example, `tau_b=tau_c`, then their fixed fields
`k(a,c)` and `k(a,b)` agree, so `a,b,c` all lie in one proper
degree-two subfield of `k(C)`, contrary to birationality.

Take any three labels in `S` and use the complementary three labels as
the normalizing points `infinity,0,1`.  Call the remaining coordinates
`a,b,c`.  The map `a:C->P^1` forgets the other two moving labels.  The
pullback of a point of `bar(M)_(0,4)` is the sum of four stable
six-point boundary divisors.  Thus

```text
[k(C):k(a)]=[k(C):k(b)]=[k(C):k(c)]=4.                (12)
```

The function `a` is fixed by `<tau_b,tau_c>`.  This automorphism group
therefore has order at most four.  Since its two generators are
distinct involutions, it has order exactly four and is a Klein four
group.  In particular `tau_b` and `tau_c` commute.  The same argument
for the other coordinates shows that the three involutions commute
pairwise.

They are also linearly independent over `F_2`.  Otherwise the group
generated by any two would be the same Klein four group; (12) would
then imply that all of `a,b,c` lie in its fixed field, again
contradicting `k(C)=k(a,b,c)`.

Consequently, any three members of `S` are linearly independent in an
elementary abelian two-subgroup of `Aut(C)`.  Such a subgroup has rank
at most three.  Indeed, after choosing an origin, its translation
subgroup lies in `C[2]` and has rank at most two, while its image in
`Aut(C,0)` has two-rank at most one.  In `F_2^3`, a set in which every
three elements are independent has size at most four: after choosing
three as a basis, the only possible fourth element is their sum.
Therefore

```text
|S| <= 4.                                               (13)
```

There is more structure in the equality case.  Because
`<tau_b,tau_c>` has order four and its fixed field is exactly the
rational field `k(a)`, it cannot consist only of translations: the
quotient of an elliptic curve by translations still has genus one.
Thus `S` contains at most one translation involution.  A four-element
cap in `F_2^3` has sum zero, so it contains an even number of
reflections.  The only possibility is that all four elements are the
reflection coset of `C[2]` inside

```text
C[2] x <one reflection>.
```

Each pair-generated Klein four group then consists of the identity,
one nonzero two-torsion translation, and two reflections.  Its rational
quotient is equivalently a two-isogeny followed by the degree-two
Kummer quotient of the resulting elliptic curve.

This equality case is itself incompatible with boundary degree one.
Choose three of its four labels as `a,b,c` and call the fourth one
`infinity`; use the other two labels as `0,1`.  Put
`tau_0=tau_infinity=tau_a tau_b tau_c`.  After forgetting `infinity`,
the configurations at `Q` and `tau_0(Q)` are projectively equivalent.
The matching projectivity fixes the labeled points `0,1`, so it has the
form

```text
M_t(x)=t x/((t-1)x+1).
```

The action of `tau_0` on `k(a)` is its nontrivial Mobius involution;
write it as `sigma_a`, and similarly for `b,c`.  The one projectivity
matching all five retained labels gives

```text
t = sigma_a(a)(1-a)/(a(1-sigma_a(a)))
  = sigma_b(b)(1-b)/(b(1-sigma_b(b)))
  = sigma_c(c)(1-c)/(c(1-sigma_c(c))).                 (14)
```

Thus `t` lies in

```text
k(a) intersect k(b) intersect k(c)=k(C)^G.
```

In particular `tau_0` fixes `t`.  Applying `tau_0` to any expression
in (14), however, interchanges the coordinate with its `sigma` image
and sends `t` to `1/t`.  Hence `t=1` or `t=-1`.  The first choice makes
all three `sigma` actions trivial and contradicts birationality.  The
second gives, for every coordinate,

```text
sigma(x)=x/(2x-1),                                     (15)
```

whose fixed values are exactly `0,1`.

The reflection `tau_0` has four fixed points, forming a torsor under
`C[2]`.  At each such point, (15) puts each of `a,b,c` in `{0,1}`.
Moreover `C[2]` acts on each coordinate either trivially or by (15),
so all four fixed points give the same triple
`(a,b,c) in {0,1}^3`.

Let `A` be the subset of `{a,b,c}` taking value zero at these points.
If `A` is nonempty, stable reduction at every fixed point contains the
boundary separating `{0} union A` from its complement.  If `A` is
empty, it contains the boundary separating `{1,a,b,c}` from
`{infinity,0}`.  These are stable boundary divisors.  The same divisor
therefore has at least four distinct points in its pullback to `C`,
contrary to degree one.  Consequently the bound improves to

```text
|S| <= 3.                                               (16)
```

The equality case in (16) has a further exact classification.  The
three independent commuting deck involutions contain at most one
translation, because the group generated by any pair has rational
fixed field.  They cannot all be reflections, as follows.

Write the three coordinates as `a,b,c`, their involutions as
`tau_a,tau_b,tau_c`, and put

```text
G=<tau_a,tau_b,tau_c>,       tau_D=tau_a tau_b tau_c.
```

If all three generators are reflections, then `tau_D` is also a
reflection and the translation subgroup of `G` is all of `C[2]`.
Let `t` generate the rational common quotient field `k(C)^G`.  Since
`k(a)`, for example, is the fixed field of
`<tau_b,tau_c>`, there is a degree-two rational map `f_a` with

```text
t=f_a(a)=f_b(b)=f_c(c).                                (17)
```

The coordinate functions and `t` extend to morphisms from the smooth
proper source curve to `P^1`; each rational map `f_i` is likewise a
morphism of projective lines. Thus (17) holds at boundary-support
points as well as on the open locus of distinct configurations.

The four fixed points of `tau_D` form a `C[2]`-torsor.  Their images in
the three coordinate lines are the ramification points
`alpha_a,alpha_b,alpha_c` of the maps in (17), and the full normalized
six-point coordinate tuple is the same at all four fixed points.  If
some `alpha` were one of `0,1,infinity`, or if two of the three
`alpha` values agreed, the same outer stable boundary cluster would
occur at all four points.  That would give its boundary pullback degree
at least four.  Hence the three `alpha` values are distinct and avoid
the anchors.

The boundary divisors where two moving labels collide with the anchor
`0` show that

```text
f_a(0)=f_b(0)=f_c(0)=A.
```

The analogous boundaries at `1,infinity` give common values `B,C`.
Send the common branch value arising from `tau_D` to infinity on the
`t`-line.  The anchor values remain finite by the preceding paragraph.
A degree-two map with double pole `alpha` and those three values is
uniquely

```text
f_alpha(x)=C+(u_alpha x+v_alpha)/(x-alpha)^2,
v_alpha=(A-C)alpha^2,
u_alpha=(B-C)(1-alpha)^2-v_alpha.                      (18)
```

The three values `A,B,C` are not all equal, since a degree-two map
cannot put the three distinct points `0,1,infinity` in one fiber.

Now use the boundary divisor separating `{a,b,c}` from the three
anchors.  At its unique pullback point, let the common limiting
coordinate be `x` and let the common value of `t` be `y`.  If `x` is
not an anchor, then `y` cannot be infinity: that would require `x` to
equal all three distinct double poles.  For finite `y`, the equation
`f_alpha(x)=y` is

```text
0=(C-y)(x-alpha)^2
  +x(B-C)(1-alpha)^2+(A-C)(1-x)alpha^2.                (19)
```

It is quadratic in `alpha`, so the three distinct `alpha` values would
force it to vanish identically.  For `x` different from `0,1`, its
linear coefficient first gives `y=B`; the constant and quadratic
coefficients then force `A=B=C`.  This is impossible.  The projective
case `x=infinity` is one of the excluded anchors.

It remains to consider a nested collision at an anchor.  At `x=0`,
the boundary separating `{a,b,c}` means that the three coordinate
expansions have the same first nonzero term: the cluster must exclude
the anchor `0`, even if it lies inside the larger cluster
`{0,a,b,c}`. Different leading terms could give that larger cluster
alone, but would not give the required inner `{a,b,c}` boundary.
Equation (17) therefore
forces the three derivatives `f'_alpha(0)` to agree.  More explicitly,
if this common leading term has order `h`, a nonzero derivative gives
order `h` for `t-A`, whereas a zero derivative gives larger order.
Thus mixed zero and nonzero derivatives are impossible; if all are
nonzero their leading coefficients must agree, and if all are zero
they already agree.  From (18),

```text
f'_alpha(0)
 =(B-A)+2(A-B)/alpha+(B-C)/alpha^2.                    (20)
```

As a polynomial in `1/alpha`, this is nonconstant unless `A=B=C`, so
one value cannot occur at three distinct `alpha` values.  The anchor
`1` is identical after replacing `x` by `1-x`.  At `infinity`, using
`1/x` as local coordinate, the derivative is

```text
(B-A)alpha^2-2(B-C)alpha+(B-C),                        (21)
```

again a nonconstant quadratic.  This excludes all nested cases and
proves that the three involutions cannot all be reflections.

Thus equality in (16) consists of exactly one translation involution
and two reflections.  The translation quotient is a smooth elliptic
anticanonical image; each reflection quotient is a singular rational
anticanonical image.  In particular, at most two anticanonical
forgetful images can be singular rational.

For a balanced degree-one six-curve, at least three forgetful images
must therefore have `delta_j=2,r_j=1`: they are birational images in
`|2L|` whose genus-one normalizations require singularity budget five.
The bound is sharp.  The quadratic-pencil construction in
[disjoint_mixed_elliptic_six_curve.md](disjoint_mixed_elliptic_six_curve.md)
has exactly three anticanonical images of the classified mixed type and
all twenty-five boundary pullbacks of degree one, supported at distinct
source points.

## 5. Scope and the next irreducible-support target

The theorem rules out the most direct extension of the five-row
construction: one cannot place a balanced fixed elliptic six-curve so
that three overlapping five-row projections are smooth anticanonical
elliptic sections.  It works for every common boundary degree `d`, not
only the minimal case `d=1`.

Balanced elliptic six-curves therefore do exist, so the geometric
genus argument cannot close the endpoint problem by itself.  Formula
(4) gives their exact projection threshold:
the normalization of `Gamma_j in |delta_jL|` has genus `g_j in {0,1}`,
while its arithmetic genus is (5), so it must carry total singularity

```text
1+5delta_j(delta_j-1)/2-g_j.                           (22)
```

For `d=1`, the only alternatives were listed after (6).  For larger
`d`, the factorization `r_jdelta_j=2d` determines the possible classes
and covering degrees.

A useful continuation would be an irreducible-support theorem showing
that a full six-row profile with a fixed genus-one support forces at
least three forgetful images to be smooth.  By (11) that would close the
route.  Equivalently, one can classify the locus where at least four of
the six projected curves acquire precisely the singularity budget
(22), and test whether an irreducible invariant relation with the full
cut support can lie there.  No such smoothness or support theorem is
currently supplied by the full-profile height calculation.

The finite combinatorial checks in
[check_six_row_elliptic_forgetful_obstruction.py](check_six_row_elliptic_forgetful_obstruction.py)
verify the twenty-five boundary partitions, the two-term pullback (2),
the divisor-class equations, the two-torsion kernel intersection, and
the degree-one involution-cap classification.
