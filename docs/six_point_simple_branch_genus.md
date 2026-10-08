# Every nondegenerate parabolic branch is simple

The degree-twelve circularity form in the
[isotropic-pencil construction](six_point_isotropic_circle_cover.md)
has a simple zero at every parameter where the six affine columns are
distinct and their unique conic is a nonsingular parabola. Consequently,
if the fixed central vector has no vanishing proper subsum, both of its
circle covers have genus exactly five. This conclusion holds for every
such central vector admitting a rational isotropic frame, rather than
only for the two previously computed examples.

This is a statement about the geometry of each fixed-weight cover. It
does not give a uniform bound for rational points or for their heights
as the weights vary.

## 1. The tangent calculation

Work over an algebraically closed field of characteristic zero. At the
parameter in question an affine change of coordinates puts the parabola
and its six points in the form

```text
y=x^2,             (x_i,y_i)=(t_i,t_i^2),
```

with the `t_i` distinct. The rowspace
`W=span(1,t,t^2)` is maximal isotropic for
`beta(v,w)=sum_i u_i v_i w_i`. Therefore

```text
S_j=sum_i u_i t_i^j=0,             j=0,1,2,3,4,
S_5 != 0.                                            (1)
```

If the last assertion failed, the invertible six-by-six Vandermonde
matrix would annihilate the nonzero vector `u`.

Choose a hyperbolic dual frame `(f0,f1,f2)` for `(1,t,t^2)`.
A local parameter `z` on the ruling through `W` is supplied by

```text
x(z)=t+z f2,          y(z)=t^2-z f1.                 (2)
```

This is the skew-matrix chart of the maximal isotropic Grassmannian.
In particular it is an actual local coordinate, not a ramified
reparametrization. At zero,

```text
beta(t^2,x')=1,       beta(t,y')=-1.                 (3)
```

Write its unique conic locally as

```text
Q_z(x,y)=c0(z)+c1(z)x+c2(z)y+A(z)x^2+B(z)xy+C(z)y^2,
Q_0=y-x^2.
```

Its coefficients can be chosen regular near zero, since the feature
matrix has rank five there; scale them so the displayed normalization
holds. Differentiate `Q_z(x_i(z),y_i(z))=0`, multiply by `u_i t_i`,
and sum. Every differentiated coefficient except `C'` multiplies one
of `S_1,...,S_4` and hence vanishes. The coordinate derivatives give

```text
0 = C'(0) S_5 - 2 beta(t^2,x') + beta(t,y')
  = C'(0) S_5 - 3.
```

Thus

```text
C'(0)=3/S_5,
(4AC-B^2)'(0)=-12/S_5 != 0.                         (4)
```

Changing the affine coordinates multiplies the quadratic discriminant
by a nonzero square constant. Changing a locally nonvanishing conic
coefficient normalization multiplies it by a local unit squared.
Neither changes the simple-zero assertion. The argument also applies
to a parameter at infinity by choosing a local chart there.

## 2. Consequence for the entire binary form

Suppose `sum_(i in S)u_i !=0` for every nonempty proper subset. The
[radius-content proof](six_point_circle_cover_radius_content.md)
then shows that the six columns are always distinct, the raw conic
cofactors have no common projective zero, and the forms `Delta=4AC-B^2`
and `K` have no common projective zero.

Every projective zero of `Delta` consequently represents a nonsingular
parabola. Its quadratic part has rank one; rank zero would make the
conic a line, contrary to six noncollinear columns and `K!=0`.
Equation (4) shows that every such zero is simple. It also rules out
`Delta` being identically zero. The raw homogeneous form has degree
twelve, so it has exactly twelve distinct zeros on the projective line
over the algebraic closure, including any zero at infinity.

The double cover `d^2=Delta(s,t)` is therefore geometrically connected
and smooth after its standard weighted-projective normalization. It is
ramified simply at these twelve points. The degree-two covering formula
`2g-2=2(-2)+12` gives `g=5`.

By [the prime-cut subsum criterion](central_subsum_cut_separation.md), this
applies eventually to a full fair six-label core with subpower
normalization and private-factor errors. It does not make its frame,
branch coordinates, or radius comparison constants independent of the
core. The original uniform arc bound remains unproved.

## 3. Verification

The exact checker `check_six_point_simple_branch.py` constructs rational
parabolic configurations, their primitive barycentric weights, and a
hyperbolic dual frame. It differentiates the conic-incidence equations
by rational interpolation and verifies `C' S_5=3` independently of
formula (4), including examples with and without vanishing proper
subsums. The general proof is the local calculation above.
