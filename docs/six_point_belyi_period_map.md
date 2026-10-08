# The integer-residue period map is Belyi

This note identifies an intrinsic rational function on each nonresonant
six-point circle cover.  It explains the degree-twelve branch divisor and
the degree-twenty conic-degeneration divisor in terms of the canonical
relative-period differential.  The resulting Belyi map has degree growing
linearly with the central integer weights, so it does not by itself give a
uniform endpoint estimate.

Fix a primitive integer vector

```text
u=(u_0,...,u_5),       sum_i u_i=0,
```

with no vanishing nonempty proper subsum.  Let `C` be either one of the two
smooth genus-five covers constructed in
[`six_point_isotropic_circle_cover.md`](six_point_isotropic_circle_cover.md):

```text
C: d^2=Delta(s,t).
```

The twelve roots of `Delta` are simple.  The degree-twenty form `K` is
squarefree, its roots are disjoint from those of `Delta`, and its ten
quadratic factors are indexed by the unordered balanced partitions
`I|I^c` of the six labels.  These facts are proved in
[`six_point_circle_cover_radius_content.md`](six_point_circle_cover_radius_content.md).

## 1. Dictionary with the isoresidual fiber

For a point of `C`, let `q_i=Z_i/Z_0` be its six normalized circle nodes.
The central moment identities are

```text
sum_i u_i q_i^j=0,             -2<=j<=2.              (1)
```

Consider on the node sphere the differential

```text
Omega_q(z)=sum_i u_i dz/(z-q_i).
```

The three nonnegative moments in (1) give a zero of order two at infinity,
and the two negative moments give a zero of order two at zero.  Equivalently,

```text
Omega_q(z)=c_q z^2 dz/product_i(z-q_i),       c_q!=0. (2)
```

Thus the fixed-`u` family is the isoresidual fiber in the stratum
`H(2,2;[-1]^6)`, after fixing the two zeros at zero and infinity and fixing
the remaining scaling by `q_0=1`.  The two isotropic rulings give its two
components.

For the source, residue convention, and genus comparison, see
[`central_moments_isoresidual_dictionary.md`](central_moments_isoresidual_dictionary.md)
and Chen--Gendron--Prado--Tahar,
[*Isoresidual curves*](https://arxiv.org/html/2412.16810v3).

## 2. The relative period is a logarithm of a rational function

Because the `u_i` are integers, the function

```text
F=product_i q_i^(u_i)                              (3)
```

is a single-valued rational function on `C`, over the same coefficient
field as the circle lift.  A primitive of `Omega_q` is
`sum_i u_i log(z-q_i)`.  Since `sum_i u_i=0`, its value at infinity is zero,
whereas its value at zero is `log F`, modulo contour periods.  Therefore the
canonical relative-period differential is exactly

```text
omega_u=-d log F.                                    (4)
```

This uses the original differential with contour residues `2 pi i u`.  If
one first rescales `Omega_q` by `1/(2 pi i)` so that the paper's residue
vector is written as `u`, then the canonical differential is instead
`-d log F/(2 pi i)`.

The period theorem says that analytic continuation changes a relative
period by an integer linear combination of the fixed residues.  It does
not say that the value of the relative period at an arbitrary point is an
integer.

There is also a direct moment check at a branch point.  Use `tau=d` as
local parameter and choose logarithms

```text
q_i=exp(r_i),       r_i=alpha_i tau+O(tau^2).
```

All `q_i` equal one at the branch.  The six `alpha_i` are distinct affine
values of a parameter on the nonsingular limiting parabola.  The five
equations in (1) say

```text
sum_i u_i exp(j r_i)=0,             j=-2,-1,0,1,2.
```

The explicit exponential interpolant

```text
P(exp(x))=(2/3)(exp(x)-exp(-x))
          -(1/12)(exp(2x)-exp(-2x))
         =x-x^5/30+O(x^7)
```

and the five moments give

```text
log F=sum_i u_i r_i
     =(tau^5/30) sum_i u_i alpha_i^5+O(tau^6).        (5)
```

The last fifth moment is nonzero.  Indeed, the `u_i` annihilate the powers
of the six distinct `alpha_i` through degree four; the Vandermonde kernel is
one-dimensional and cannot also annihilate degree five.  Thus `F-1` has
exact order five and `d log F` has exact order four at every branch point.

## 3. Exact divisor of the canonical differential

The homogeneous differential

```text
eta=d^3 (s dt-t ds)/K                              (6)
```

is well-defined on `C`: its numerator and denominator have the same
homogeneous weight twenty.  At a simple root of `Delta`, take `d` as local
parameter.  The base parameter then has order two, so (6) has order
`3+1=4`.  At each of the two lifts of a simple root of `K`, it has a simple
pole.  In any projective chart it is regular and nonzero away from the
`Delta`- and `K`-divisors, so there is no additional contribution from the
chosen parameter infinity.  Hence on either genus-five component its divisor
consists of

```text
12 zeros of order 4,        40 simple poles.          (7)
```

This has degree `48-40=8=2 genus(C)-2`.

The proportionality with (4) can be computed exactly.  Work in an affine
parameter `r` on the pencil and put

```text
a_i=2A x_i+B y_i+c1,
b_i=Delta y_i+delta,       delta=2A c2-B c1,
Z_i=d a_i+i b_i,           Z_i_tilde=d a_i-i b_i.
```

Then `Z_i Z_i_tilde=4AK`.  Since `sum_i u_i=0`, normalization by `Z_0`
cancels from the logarithmic derivative, and

```text
d log F/dr=sum_i u_i Z_i'/Z_i.
```

Write `beta(v,w)=sum_i u_i v_i w_i`.  Expanding
`sum_i u_i Z_i' Z_i_tilde`, the real terms vanish by
`beta(a,a)=beta(b,b)=0`.  The term containing `d'` vanishes by
`beta(a,b)=0`, and differentiating that equality gives
`beta(a',b)=-beta(a,b')`.  Hence

```text
d log F/dr = i d beta(a,b')/(2AK)
           = i d^3 beta(x,y')/K.                     (8)
```

The last equality uses
`beta(a,b')=Delta beta(a,y')=2A Delta beta(x,y')`.
For a homogeneous linear pencil, `beta(x,dy)` is a fixed nonzero multiple
`kappa(s dt-t ds)`.  In the standard chart
`x=e1+r f2`, `y=e2-r f1`, one has `kappa=-1`.  Therefore

```text
d log F = i kappa d^3(s dt-t ds)/K.                  (9)
```

This proves the differential identity directly, including its constant.

## 4. Belyi values and exact degree

At a root of `Delta`, formula (3) has `q_i=1` for all six labels.  Therefore
`F=1`, and the order-five conclusion above says that the twelve branch
points have ramification index five over `1`.

A root of the quadratic factor indexed by `I|I^c` is the degeneration into
two lines carrying the two triples.  Choose an affine chart with `A!=0` at
this root.  The two nonparallel component lines have equations `Z=0` and
`Z_tilde=0`, and none of the six marked points is their intersection.  On
one lift, `Z_i` vanishes for the three labels of one component and is a unit
for the other three.  Since `Z_i Z_i_tilde=4AK` and the root of `K` is
simple, each of those three zeros is simple.  Normalization by `Z_0`
cancels because `sum_i u_i=0`, so

```text
ord(F)=plus-or-minus sum_(i in I)u_i.                (10)
```

The order is nonzero by the nonresonance hypothesis.  The two choices of
`d` are related by circle inversion, which sends `F` to `F^(-1)` and
reverses the sign in (10).  Thus the two lifts of each base root contribute
one zero and one pole, each of order `|sum_(i in I)u_i|`.

There are no other zeros or poles: at either one, `d log F` would have a
simple pole with nonzero integer residue, whereas (9) has poles only over
`K=0`.

There are two roots for each of the ten balanced partitions.  It follows
that

```text
degree(F)=2 sum_(unordered balanced I|I^c)
                    |sum_(i in I)u_i|.               (11)
```

Away from the zeros and poles of `F`, its critical points are precisely the
zeros of `d log F`.  Equations (6)--(10) show that all of them lie over `1`.
Consequently

```text
F:C -> P^1 has critical values contained in {0,1,infinity}.
```

It is a Belyi map.  Formula (11) also states its limitation for the endpoint
problem: its degree is unbounded and grows linearly in the proper subset
sums of the moving central vector.  The integrality of its monodromy
periods does not make the circle-node positions, angular gaps, or branch
parameter values integral.

## Appendix. Rational descent and its exact ramification data

Use an affine coordinate `r` on the rational isotropic ruling, so the cover
has equation `d^2=Delta(r)` with coefficients in `Q`. Put

```text
G=(F+F^(-1))/2,
H=(F-F^(-1))/(2i d).                                (12)
```

Both `G` and `H` belong to `Q(r)`. Indeed the hyperelliptic involution
`sigma:d->-d` sends each normalized node to its inverse and hence sends
`F` to `F^(-1)`. It fixes both expressions in (12): the numerator and
denominator of `H` both change sign. Thus `G,H` first descend to `Q(i)(r)`.
The coefficient automorphism `i->-i`, fixing `r,d`, also sends `F` to
`F^(-1)`, because `Z_i_tilde=4AK/Z_i`. It again fixes both expressions,
proving their descent to `Q(r)`.

Their exact identity over `Q` is

```text
G^2+Delta H^2=1,       F=G+i d H.                    (13)
```

In particular, over `Q(i)` the rational function `G^2-1` has squareclass
`Delta`. The factor of `i` in (12) is essential for the plus sign in (13).
Neither `H` nor `F-F^(-1)` is identically zero, since `F` is nonconstant.

Write `D=degree(F)`, given by (11). The base rational map `G:P^1_r->P^1`
has degree `D`: its pullback to `C` has degree twice its base degree, while
`(F+F^(-1))/2` has degree `2D` as the composition of `F` with a degree-two
rational map. Its fibers over the three possible critical values are:

| Value of `G` | Points on the base | Ramification indices |
|---|---:|---|
| `+1` | twelve roots of `Delta` | `5` at each |
| `+1` | `(D-60)/2` additional points | `2` at each |
| `-1` | `D/2` points | `2` at each |
| infinity | twenty roots of `K` | `|sum_(i in I)u_i|` at a root belonging to `I|I^c` |

Each balanced partition contributes two roots in the last row. Here and
in the table all roots are counted geometrically over the algebraic closure.

To verify the finite fibers, use

```text
G-1=(F-1)^2/(2F),       G+1=(F+1)^2/(2F).            (14)
```

At a branch point, `F-1` has order five on `C`, while the base projection
has local degree two. Hence `G-1` has order five on the base. Every other
point with `F=1` is unramified for `F` and for the base projection, and
comes with a distinct partner under `sigma`; (14) gives a base point of
index two. The degree of the `F=1` fiber then gives `(D-60)/2`. The same
argument at `F=-1`, where there are no branch points, gives `D/2` points of
index two. In particular these counts imply `D>=60`.

At a root of `K`, the two points of `C` above it are respectively a zero
and a pole of `F` of the same order. Formula (12) therefore gives a pole
of that order for `G` on the base. There are no other poles. The composition
description, or (14) together with the known critical points of `F`, shows
that no other critical values occur. After an affine change of target,
`G` is therefore a rational Belyi map defined over `Q`.

As an independent ramification check, the contributions are

```text
12*4+(D-60)/2+D/2+(D-20)=2D-2,
```

exactly Riemann--Hurwitz for a degree-`D` map of projective lines. The
ramification data do not identify a particular dessin or prove that its
underlying graph is simple. The identity (13) and the fixed number of
index-five points provide additional structure, but the degree and the
pole orders still depend on the moving central weights; no new uniform
arithmetic separation is asserted here.
