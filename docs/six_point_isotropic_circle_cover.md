# An exact fixed-weight circle cover of the isotropic pencil

For a fixed nonzero rational six-vector `u` with `sum u_i=0` and hyperbolic
`diag(u)`, the circle configurations with central vector proportional to `u`
admit an explicit description on two isotropic pencils. Each pencil with a
nonempty nondegenerate circle locus is lifted to a double cover of its
projective line, of genus at most five. Two exact examples have squarefree
degree-twelve branch polynomials, hence genus five; the cover is not
generically an elliptic quartic.

The subsequent [simple-branch proof](six_point_simple_branch_genus.md)
shows more: both covers have genus exactly five whenever `u` has no
vanishing proper subsum. The [prime-cut criterion](central_subsum_cut_separation.md)
verifies that condition on an eventual full fair six-label core.

The construction fixes the vector `1`, retains the affine chart, and gives
an exact Gaussian lift with its common gcd. Thus no column rescaling or
projective denominator is silently discarded. It supplies a parametrization
and an upper height estimate, not a lower radius bound excluding a fair
six-label core.

The later [exact radius-content calculation](six_point_circle_cover_radius_content.md)
replaces the crude upper estimate below by an exact degree-twenty expression
with all contents specified. Under its fixed-weight nonvanishing-subsum
hypothesis it proves `N asymp H^20` and reduces bounded endpoint scale to
`d^2=Delta(s,t)` with `|d|=O(H)`. Its comparison constants still depend on the
fixed frame.

## 1. Two rational pencils preserving the zeroth moment

Write `beta(v,w)=sum_i u_i v_i w_i`, and let `e0=1`. Choose a rational
hyperbolic frame

```text
(e0,e1,e2; f0,f1,f2),
beta(e_i,e_j)=beta(f_i,f_j)=0, beta(e_i,f_j)=delta_ij.
```

Such a frame exists because the form is split and `e0` is nonzero isotropic.
The maximal isotropic three-spaces containing `e0` correspond to maximal
isotropic two-spaces in `e0^perp/<e0>`. Their two rulings are

```text
W_plus(s:t)  = span(e0, s e1+t f2, s e2-t f1),
W_minus(s:t) = span(e0, s e1+t e2, s f2-t f1).        (1)
```

All three generators are independent for every `(s:t)` and are pairwise
orthogonal. These are complete pencils: a plane transverse
to `span(f1,f2)` in the split four-space is the graph of a skew two-by-two
matrix, giving the first ruling, and exchanging `e2,f2` gives the second.
The remaining member of each ruling is its parameter at infinity, `s=0`.

One can construct the dual frame without hiding denominators. Given the
three-by-six row matrix `E=(e0,e1,e2)`, choose any rational `F_temp` with
`E diag(u) F_temp^T=I`. Put

```text
B0=F_temp diag(u) F_temp^T,
F=F_temp-(1/2) B0 E.                                (2)
```

Then `E diag(u) F^T=I` and `F diag(u) F^T=0`. Every rational entry created
by (2) belongs to the frame height used below.

## 2. The conic and its degree-twelve circularity cover

For either pencil, write the last two row vectors as `x(s,t),y(s,t)`.
The six affine points `(x_i,y_i)` have Veronese feature rows

```text
F_i=(1,x_i,y_i,x_i^2,x_i y_i,y_i^2).
```

Because the three rows of (1) are totally isotropic,

```text
sum_i u_i F_i = 0.                                 (3)
```

Thus the six feature rows have rank at most five. On the open locus of six
distinct columns and feature rank five, their unique conic is

```text
c0+c1 x+c2 y+A x^2+B xy+C y^2=0.                    (4)
```

Take signed five-by-five cofactors of any five independent feature rows.
The coefficient degrees in the homogeneous pencil parameters `(s,t)` are

```text
deg(c0,c1,c2,A,B,C) = at most (8,7,7,6,6,6).         (5)
```

These degree counts include every denominator after one common clearing of
the frame. They follow by adding the column degrees `0,1,1,2,2,2` and omitting
the cofactor's column. Removing a common polynomial factor is permissible,
but must be done explicitly; the examples below have none.

There is a useful check on possible common zeros. If two affine columns
`i,j` coincide, the vector `a=e_i/u_i-e_j/u_j` (here `e_i` denotes a coordinate
unit vector) is orthogonal to `W`. A maximal isotropic space is its own
orthogonal complement, so `a` lies in `W` and
`beta(a,a)=1/u_i+1/u_j=0`. Thus a collision requires `u_i+u_j=0`.
If no two weights are opposite, all columns stay distinct throughout both
pencils, even over the algebraic closure. Moreover four columns cannot be
collinear: a linear form vanishing on those four belongs to `W`, and its
orthogonality to all of `W` would make the other two projective columns
proportional. Hence every five columns determine a unique conic, and the
cofactor vector has no common zero on the parameter line. This verifies
that removing an unproved common factor is particularly unjustified for
weights without opposite pairs.

Set

```text
Delta(s,t)=4A(s,t)C(s,t)-B(s,t)^2.
```

The required cover is

```text
d^2=Delta(s,t),                                    (6)
```

where `d` has homogeneous weight six when the raw cofactors are used.
Its branch polynomial has degree at most twelve.
If `Delta` is identically zero, that pencil has no `d!=0` circle locus and
is discarded. Otherwise square factors are removed
in forming the normalization; every resulting geometrically integral double
cover has genus at most five. If `Delta` is a square, the cover splits into
rational components. If the degree-twelve binary polynomial is squarefree,
there are twelve branch points and Riemann--Hurwitz gives genus five.

The condition in (6) retains exactly the two nonreal directions at infinity:
the quadratic part of (4) has discriminant `B^2-4AC=-d^2`. Merely knowing that
the conic is rational does not impose this condition. All degenerate conics,
coincident columns, and `d=0` are excluded from the circle locus.

## 3. The exact circle map preserves the weights

Suppose `(s:t:d)` is rational, `d!=0`, and (4) is nonsingular with the six
distinct rational points just constructed. Then `Delta=d^2>0`, so the
quadratic part is definite, with `A!=0`. Define rational Gaussian numbers

```text
Z_i = d(2A x_i+B y_i+c1)
      + i(d^2 y_i-B c1+2A c2).                      (7)
```

Direct substitution using (4) gives the same positive norm for every row:

```text
M = |Z_i|^2
  = 4A(C c1^2-B c1 c2+A c2^2-c0 d^2).              (8)
```

For a geometric verification, the conic center is obtained from
`2A xc+B yc=-c1`, `B xc+2C yc=-c2`. Writing `xi=x-xc`, `eta=y-yc`,
the linear map

```text
(xi,eta) -> (2A xi+B eta)+i d eta
```

has norm `4A(A xi^2+B xi eta+C eta^2)`. Formula (7) multiplies this map
by `d` and clears the center's denominator. Since the conic has real points
and is nonsingular, its common norm is positive.

The rows `1,Re(Z),Im(Z)` are a rational change of basis of the same space
`W`. Therefore they retain the exact central identities

```text
sum_i u_i Z_i^j=0,      j=0,1,2.                    (9)
```

There is no separate scaling of the six columns. The unit-circle nodes are
`q_i=Z_i/Z_0`. The two choices of `d` give conjugate normalized configurations.

Conversely, any six distinct rational unit-circle nodes with central vector
`u` give a space `span(1,Re(q),Im(q))` in one of the pencils (1). Changing
between its affine bases sends the unit circle to (4); the quadratic part's
determinant is a rational square up to the square of the conic's common
coefficient scalar. Thus (6) has a rational lift, and (7) recovers the
original circle configuration up to common rotation and reflection.

## 4. Least radius and an explicit height cost

Take one positive integer `D` clearing all coordinates of (7), and let
`G=gcd_G(D Z_0,...,D Z_5)`. The primitive tuple and squared radius are exactly

```text
z_i=D Z_i/G,
N=R^2=D^2 M/N_G(G).                                (10)
```

This is independent of the chosen clearing denominator. A Bezout combination
of the primitive coordinates shows that any other Gaussian integer realization
of this same Gaussian rational projective tuple is an integral Gaussian
multiple. Thus (10) is its least radius, not just a convenient upper bound.
It retains all singleton, pair, and balanced conductor contributions, even
though only the pair squareclasses appear in the Witt test.

An independent formula uses every pair. Put

```text
c_ij=Im(Z_i conjugate(Z_j))/(M-Re(Z_i conjugate(Z_j))),
L=lcm of the denominators of all c_ij, Q_ij=L c_ij.
```

For distinct rows the displayed denominator is positive; antipodal rows give
`c_ij=0`. Set `g_ij=gcd(Q_ij,L)`, and `epsilon_ij=2` when both reduced
entries `Q_ij/g_ij,L/g_ij` are odd, otherwise one. Then

```text
N=lcm_(i<j) (Q_ij^2+L^2)/(g_ij^2 epsilon_ij),        (11)
```

by the [all-edge radius theorem](integer_cotangent_lcm_height_target.md).

For an explicit height audit, clear the frame first so that the coefficients
of all linear `x_i(s,t),y_i(s,t)` are integers of absolute value at most `T`.
At coprime integer parameters let `H=max(|s|,|t|)` and `K=2TH>=1`.
Raw cofactor bounds give

```text
|A|,|B|,|C| <=120 K^6,
|c1|,|c2| <=120 K^7, |c0|<=120 K^8.
```

At a rational point of (6), the raw integer `Delta` is a rational square,
so `d` is an integer. Consequently (7) is integral already, and the crude
fully explicit bounds

```text
max_i max(|Re Z_i|,|Im Z_i|) <=10^6 K^13,
N <=2*10^12 K^26                                 (12)
```

hold before or after primitive division. The factor `T` includes the actual
height of the chosen isotropic frame and all its cleared denominators.
Neither `T` nor `H` is controlled by the pair squareclasses alone. Equation
(12) is an upper bound and supplies no positive radius saving on a full
fair profile. In particular it does not justify treating a curve with
central coefficient height `exp(7w+o(w))` as a fixed curve of subpower height.

## 5. Two exact genus-five witnesses

Starting from the rational circle parameters `0,1,2,3,4,5` gives the primitive
central vector

```text
(1,-20,250,-1000,1445,-676).
```

The exact frame construction (2), first pencil (1), and cofactor construction
give degrees `(8,7,7,6,6,6)`, common polynomial gcd one, and a degree-twelve
polynomial `Delta` with `gcd(Delta,Delta')=1`. The parameters
`0,1,2,4,7,11` independently give the same degree and squarefreeness conclusions
with central vector

```text
(-45,616,-3850,15895,-27500,14884).
```

For each of these two first-ruling witnesses, `Delta(t)` is irreducible
of degree twelve over `Q` and has exactly two real roots. The SymPy checker
`check_six_point_branch_arithmetic.py` certifies these statements with exact
`factor_list` and `Poly.count_roots(-oo,oo)` calls. Thus the branch roots need
not all be real. These two examples do not establish a general real-root
count for other weights or rulings.

Run `python3 docs/check_six_point_isotropic_circle_cover.py`. The standard-
library checker constructs all coefficients over `Q`, verifies the polynomial
isotropy and conic equations, performs exact polynomial gcds, and checks
(7)--(11) on the rational point coming from each initial configuration.
The degree-twelve squarefreeness certifies genus five without a numerical
root or genericity assumption. Special weights, including the earlier
[fixed-weight elliptic family](six_point_fixed_moment_elliptic_obstruction.md),
can give exceptional covers and require their own factor audit.
