# A quadratic-cover reduction for the ten-factor moment code

The ten-factor code of the
[certified real outer-cut witness](tenfactor_real_outercut_witness.md)
has two identical columns. Swapping their coefficients is an involution.
On the nondegenerate coefficient locus, its quotient has an explicit
description by three homogeneous cubics in five variables. Recovering
rational coefficients from a rational quotient point requires one
specified discriminant to be a rational square.

This is an exact inverse correspondence, not merely a necessary
elimination condition. It supplies neither a rational point on the
nondegenerate locus nor a global component or genus computation.

## The linear kernel and its involution

Use the sign matrix `S` and `A=(S_1-S_0,...,S_4-S_0)/2` of the witness,
with zero-based indices. Its coefficient columns four and seven are
identical. Solving `Aa=0` gives exactly

```text
a_0=a_4-a_5+a_6+a_7,
a_1=a_4-a_5+a_7-a_8+a_9,
a_2=a_5-a_6+a_8-a_9,
a_3=2a_4-a_5+2a_7-a_8,
```

with `a_4,...,a_9` free. Set

```text
t=a_4+a_7,   b=a_5,   c=a_6,   d=a_8,   e=a_9,
x=a_4,      a_7=t-x,
u=t-b+c,    v=t-b-d+e,    w=b-c+d-e,    q=2t-b-d.
```

Then the full coefficient vector is

```text
a=(u,v,w,q,x,b,c,t-x,d,e).                             (1)
```

The involution is `x -> t-x`, keeping all five quotient coordinates
fixed. Distinct absolute values imply `a_4!=-a_7`, hence `t!=0`.

## Three cubic equations and the exact inverse

Define

```text
H_1=u^3-v^3-c^3-d^3+e^3,
H_2=u^3+v^3+w^3-q^3-d^3,
H_3=b^3-c^3+d^3-e^3-w^3,
B=v^3-q^3+t^3-e^3.                                    (2)
```

On `t!=0`, put

```text
p=B/(3t),                D=t^2-4p.                   (3)
```

The complete original odd-moment system `Aa=0,Aa^3=0` is equivalent to

```text
H_1=H_2=H_3=0,             x^2-tx+p=0.               (4)
```

To verify both directions, the linear equations vanish identically in
(1). If `F_0,...,F_3` are the four entries of `Aa^3`, direct expansion
gives the polynomial identities

```text
F_0=B+3t(x^2-tx),
F_1=F_0+H_1,
F_2=H_2,
F_3=H_2+H_3.                                          (5)
```

Thus (4) is sufficient as well as necessary. All three `H` polynomials
are homogeneous cubics; `p` is homogeneous of degree two. Projectively,
this describes a quadratic cover of the open part `t!=0` of their common
zero locus in `P^4`.

For rational quotient coordinates, a rational lift exists precisely when
`D` is a square in `Q`; the two choices are

```text
x=(t+sqrt(D))/2,             t-x=(t-sqrt(D))/2.        (6)
```

To obtain the required nondegenerate moment solution one must additionally
check that **all ten** entries in (1) are nonzero with pairwise distinct
absolute values. In particular `p!=0` and `D!=0` are necessary, but they
do not replace the other zero/collision exclusions. A real nondegenerate
lift requires `D>0`; positivity of the ten magnitudes then follows by
absorbing each coefficient sign into its source column.

The real witness therefore does not imply rational flexibility. It proves
existence of a suitable real quotient point and a positive real
discriminant, not a rational quotient point with square discriminant.

## Lower-dimensional arithmetic structure

The third equation does not involve `t`. It is the first-and-third-moment
equation for the five signed values

```text
(b,-c,d,-e,-w),              b-c+d-e-w=0.             (7)
```

The first equation is quadratic in `t`; its leading coefficient is
`3(c+d-e)=3(u-v)`, which is nonzero on the admissible locus because
`a_0!=a_1`. There are useful compact factorizations. Put
`h=c+d-e`, `k=c-d+e`, and `y=t-b`. Then exactly

```text
H_1/3=h y(y+k)-de(d-e),
H_3/3=k b(b-k)-(c+e)(d-e)(c-d).                       (8)
```

Here `k=b-w` is also nonzero on the admissible locus. In particular
rational quotient points must make both quantities

```text
D_y=k^2+4de(d-e)/h,
D_b=k^2+4(c+e)(d-e)(c-d)/k                            (9)
```

rational squares. These two square conditions still leave `H_2=0` and
the lifting square condition `D=t^2-4p`; they are not sufficient alone.
Thus one can first study the explicit cubic surface (7),
then impose the compatibility of that quadratic with `H_2`, and finally
the square condition (3). These reductions do not assert that the
resulting curve is rational, irreducible, or smooth globally. In
particular, elimination can introduce components on excluded coordinate
or collision loci, and a generic complete-intersection genus is not
a genus computation for the admissible component.

There is a rigorous local smoothness statement at the certified witness.
Its four-by-four cubic Jacobian in the fixed chart `a_5=1,a_9=-2/5` is
invertible by the contraction certificate. In that box,
`t>0.75` and `|2x-t|>0.85`. Hence
`partial F_0/partial x=3t(2x-t)` is nonzero. Locally one can eliminate
`x` from (5), and the three quotient equations have rank three. The
quotient is consequently a smooth curve near the witness after
projectivization, and the quadratic cover is unramified there. This is
local information only.

The [checker](check_tenfactor_involution_quotient.py) verifies the full
linear kernel and all four identities (5) with sparse rational
polynomials, (8), the homogeneous degrees, the quadratic leading coefficient,
and the explicit nonvanishing margins in the certified box.

## A bounded rational-point diagnostic

The separate [bounded search](check_tenfactor_quotient_bounded_search.py)
enumerates primitive integer representatives of the critical quintet
`(b,-c,d,-e,-w)`. The
[critical five-root sign theorem](critical_odd_moment_root_order.md)
puts their ordered
absolute values `r_1<...<r_5` in the signed pattern `++--+`, up to global
sign. If `s=r_1+r_2`, their two moment equations imply

```text
r_3+r_4=s+r_5,
r_3 r_4=s r_5+r_1 r_2 s/(s+r_5).
```

Integer divisibility and a square discriminant therefore enumerate all
such quintets at a stated bound. The search tries all 120 assignments to
the five signed roles, all rational roots of `H_1` in `t`, and then
tests `H_2` exactly. It allows unbounded rational values for the remaining
coefficients once the five selected roles have been fixed.

At maximum primitive quintet magnitude 300, there are 593 quintets and
232 assignments with rational quadratic roots; none satisfies `H_2`.
At bound 1000 the corresponding counts are 3645 and 444, again with no
`H_2` solution. These are finite diagnostics for this particular code and
this normalized five-coordinate height. They do not prove that the curve
has no nondegenerate rational points, and they supply no obstruction for
arbitrary circle tuples or other ten-factor codes.
