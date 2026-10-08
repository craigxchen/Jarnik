# An explicit balanced elliptic six-curve

The remaining degree-one case in
[six_row_elliptic_forgetful_obstruction.md](six_row_elliptic_forgetful_obstruction.md)
actually occurs.  This note gives an explicit curve over `Q` with all
twenty-five boundary pullbacks of degree one.  Three forgetful images
are anticanonical: one is smooth elliptic and two are singular rational.
Thus the fixed-genus geometric route cannot exclude every balanced
six-point elliptic support.  The twenty-five incidences occur at only
nineteen source points: six points lie on two disjoint boundary
components.  This overlap means the example is not yet a full
independent-prime profile.

## 1. Three maps in one quadratic pencil

Put

```text
f_lambda(x)=(1+2lambda)x(x-1)/(1+lambda x)              (1)
```

and choose

```text
lambda_a=3,       lambda_b=-57/121,       lambda_c=-9/25. (2)
```

For two parameters, the homogeneous cross-product is

```text
P_lambda Q_mu-P_mu Q_lambda
 =-(lambda-mu)X(X-Y)Y(X-2Y).                            (3)
```

Consequently all three maps take common values at the four common
arguments

```text
x=0,1,infinity,2.                                      (4)
```

Their finite critical points satisfy

```text
lambda x^2+2x-1=0.                                     (5)
```

For the three parameters in (2), the critical points and values are

| map | critical points | branch values |
|---|---|---|
| `f_a` | `1/3,-1` | `u=-7/9, v=-7` |
| `f_b` | `11/3,11/19` | `u=-7/9, s=-7/361` |
| `f_c` | `5,5/9` | `v=-7, t=-7/81` |

The four branch values are distinct.  The branch pairs therefore have
the required path form

```text
f_a:{u,v},       f_b:{u,s},       f_c:{v,t}.            (6)
```

## 2. The genus-one source and its involutions

Let `E` be the normalization of

```text
f_a(a)=f_b(b)=f_c(c)=T.                                (7)
```

Equivalently, each coordinate satisfies

```text
(1+2lambda_i)x_i^2
 -((1+2lambda_i)+lambda_i T)x_i-T=0.                   (8)
```

The three discriminant squareclasses have branch supports
`{u,v}`, `{u,s}`, and `{v,t}`.  They are independent: the branch value
`s` occurs only in the second, `t` only in the third, and then the first
remains.  Hence (7) is a connected degree-eight
`(Z/2Z)^3`-cover of the `T`-line.

The inertia sign vectors at `u,v,s,t`, in the order `a,b,c`, are

```text
(1,1,0),       (1,0,1),       (0,1,0),       (0,0,1).   (9)
```

Each branch value has four points of ramification index two.  Thus

```text
2g(E)-2=8(-2)+4*4=0,
```

so `E` has genus one.  The fiber over `T=0` is unramified and consists
of the eight rational choices `a,b,c in {0,1}`.  In particular `E` has
a rational point and is an elliptic curve over `Q` after choosing one.

The involution changing only the `a` root never occurs as an inertia
vector in (9), so it is fixed-point-free and is a translation.  The
individual `b` and `c` involutions do occur and are reflections.  Thus
the three degree-two forgetful maps have exactly the mixed type left by
the abstract classification.

The functions `a,b,c` generate the function field of `E`: each recovers
`T` by (1), and each quadratic square root in (8) is linear in its
coordinate.  The resulting map

```text
E -> bar(M)_(0,6),       Q |-> (infinity,0,1,a(Q),b(Q),c(Q)) (10)
```

is therefore birational onto its image.

## 3. The twenty-five boundary contacts

All contacts occur over `T=0, infinity, 2`.

At `T=0`, every coordinate independently equals `0` or `1`.  For each
nonempty subset `S` of `{a,b,c}` taking value zero, there is the outer
cluster `{0} union S`; these give seven boundary divisors.  The seven
nonempty subsets taking value one similarly give clusters
`{1} union S`.  A mixed binary pattern lies on the two disjoint outer
boundary components, one at zero and one at one.  The eight points in
this fiber therefore give

```text
7+7=14                                                     (11)
```

distinct boundary contacts.

At `T=infinity`, the two inverse values of each map tend to `infinity`
and to its other pole.  The three other poles are

```text
-1/3,       121/57,       25/9,                         (12)
```

which are distinct nonanchors.  Every nonempty subset of the three
coordinates may choose the infinite root, giving the seven clusters
`{infinity} union S`.

Finally, all three maps take the value `T=2` at `x=2`.  Their other
inverse values are

```text
-1/7,       -121/7,       -25/7.                        (13)
```

They are distinct nonanchors.  A subset of coordinates collides at
`x=2` precisely when it has size at least two.  This gives the three
moving pair boundaries and the one moving triple boundary, hence four
contacts.  In total,

```text
14+7+4=25.                                              (14)
```

These are exactly the twenty-five stable bipartitions, each once.  They
are supported at nineteen points of `E`: all eight points over `T=0`,
seven of the eight points over `T=infinity`, and four of the eight
points over `T=2`.  Each of the six mixed binary patterns over `T=0`
lies on two disjoint boundary components, accounting for the six extra
incidences.
There are no omitted collisions: (3) says that two moving coordinates
can agree only at the four arguments in (4), and a moving coordinate
can meet an anchor only in the already listed fibers.

## 4. Contact orders

The three special target values `0,2,infinity` are not branch values
in (6), so they are local parameters on every point of their fibers.
At the common roots `0,1,2`, the inverse derivatives for the three
maps are respectively

```text
x=0:       -1/7, -121/7, -25/7,
x=1:        4/7,   64/7,  16/7,
x=2:        7/15,  7/135, 7/39.                        (15)
```

At infinity, using `1/x` and `1/T` as local coordinates, the three
slopes are

```text
7/3,       -7/57,       -7/9.                          (16)
```

Every displayed list has distinct nonzero entries.  Thus each outer
cluster has smoothing parameter of order one, its marked points are
distinct on the first bubble, and no unlisted nested boundary occurs.
Every boundary pullback in (14) consequently has degree exactly one.

This is equality of individual geometric boundary degrees.  A full
independent arithmetic profile requires different cut cores to be
pairwise coprime after controlled trimming.  Boundary functions meeting
at the same point of `E` share a local uniformizer, so the six double
incidences create common local valuation support along an elliptic
orbit.  Nothing here removes that support at subquadratic cost.

## 5. A non-torsion rational orbit exists

The three completely split unramified fibers over `T=0,2,infinity`
give twenty-four distinct rational points of `E`.  One can certify that
some are non-torsion without invoking a classification of rational
torsion.

In the sign-vector notation of (9), the kernel of the character
`chi_b chi_c` is

```text
{1,tau_a,tau_b tau_c,tau_a tau_b tau_c}.
```

All three nonidentity elements are translations: `tau_a` is free, the
product of the two distinct commuting reflections is a nonzero
translation, and their product is the third one.  This kernel is
therefore the full translation subgroup `E[2]`.  The corresponding
quadratic quotient has the model

```text
Y^2=(T+7)(T+7/9)(T+7/81)(T+7/361).                    (17)
```

After choosing any rational point as origin, quotienting by this
kernel is the multiplication-by-two quotient and is `Q`-isomorphic to
`E`.  The good reductions of (17) have respectively 20 and 24 points
over `F_13` and `F_17`.  Rational torsion injects prime-to-13 into the
first reduction and prime-to-17 into the second; using the other prime
for the possible 13- or 17-primary part shows that its order divides
`gcd(20,24)=4`.  Since `E(Q)` already contains twenty-four distinct
points, it has positive rank.  Choosing a non-torsion one gives an
actual elliptic orbit on this balanced six-curve.

The exact arithmetic and the complete bipartition enumeration are
checked by
[check_mixed_elliptic_six_curve_example.py](check_mixed_elliptic_six_curve_example.py).
