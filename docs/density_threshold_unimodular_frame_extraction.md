# The `R^(1/6)` threshold for bounded unimodular frame extraction

This note isolates what short-subarc averaging can and cannot supply for an
endpoint cluster.  From sufficiently many points one can extract five
consecutive points of bounded affine type, in a genuine unimodular coordinate
frame.  The required density is of order `R^(1/6)`.  Thus this mechanism gives
a density-sensitive bound, but it cannot turn the mere failure of a uniform
point-count bound into a fixed quadrilateral.

## 1. Statement

Fix `C>0`.  Let `P_0,...,P_(N-1)` be distinct integer points, in order, on an
arc of a circle centered at the origin and of radius `R`.  Suppose the arc
length is at most `C sqrt(R)`.

For every `K>=1`, if

```text
N-4 >= (32 C^3/K)^(1/3) R^(1/6),                    (1)
```

then five consecutive points admit a description

```text
P_(j+r)=P_j+T u_r,       0<=r<=4,                    (2)
```

where `T in GL_2(Z)`, `u_0=0`, and all `u_r` belong to a finite set depending
only on `K`.  In particular their coordinate heights are bounded solely in
terms of `K`.

Combining this extraction with
`fixed_quadrilateral_unimodular_isolation.md` gives a threshold `R_0(C)` such
that

```text
N < 4+32^(1/3) C R^(1/6).                            (3)
```

for `R>=R_0(C)`. This is much weaker than the already established
`(1+epsilon) log R/log log R` estimate. Its purpose is only to identify the
power loss in elementary averaging; it is not a growth-rate improvement.

## 2. A five-point window with small determinants

Let `alpha_i>0`, `0<=i<=N-2`, be the angular gaps.  Their sum `Delta` obeys

```text
Delta <= C/sqrt(R).                                   (4)
```

There are `N-4` windows consisting of four consecutive gaps.  The sum of all
window sums is at most `4 Delta`, since each gap occurs in at most four
windows.  Hence some window has total angular width

```text
S <= 4 Delta/(N-4).                                   (5)
```

Write its five points as `Q_0,...,Q_4`, its four chord vectors as
`e_i=Q_(i+1)-Q_i`, and

```text
D_i=|det(e_(i-1),e_i)|,       1<=i<=3.
```

For two successive angular gaps `alpha,beta`, the exact circle formula is

```text
D=4R^2 sin(alpha/2) sin(beta/2) sin((alpha+beta)/2).
```

Using `sin x<=x` and `alpha+beta<=S` gives

```text
D <= (R^2/2) alpha beta(alpha+beta) <= (R^2/2)S^3
  <= 32 C^3 sqrt(R)/(N-4)^3.                          (6)
```

No three distinct points of a circle are collinear, so every `D_i` is a
positive integer.  Under (1), all three lie in `{1,...,K}`.

## 3. Bounded determinants give a genuine unimodular type

For completeness, the bounded-determinant conclusion can be upgraded from a
rational frame to (2).  The argument in
`five_point_bounded_triangle_determinants.md`, through its equation (5), shows
that every determinant `det(e_i,e_j)` is bounded by

```text
B(K)=2^18 K^10.                                       (7)
```

Set `A=[e_0 e_1]` and `d=|det A|`.  Then `1<=d<=K`.  Left multiplication by
some matrix in `GL_2(Z)` puts `A` into column Hermite normal form.  There are
only finitely many such forms with determinant at most `K`; in particular all
their entries are bounded in terms of `K`.

For each partial sum `p_r=Q_r-Q_0`, the vector

```text
U_r=d A^(-1)p_r
```

is integral.  Its coordinates are determinants of `p_r` with `e_0,e_1`, hence
are bounded using (7) and the fact that `p_r` is a sum of at most four chords.
After the same unimodular left multiplication that reduces `A`,

```text
p'_r=A' U_r/d.
```

Both `A'` and `U_r` range over finite sets, so the integral vectors `p'_r`
range over a finite set depending only on `K`.  Taking the inverse of that
left multiplication as `T` proves (2).  The first four nodes have no three
collinear, as they lie on a circle.

## 4. Consequence of fixed-type isolation

Take `K=1`.  If (1) holds, Section 3 places the extracted five-point block in
one of finitely many unimodular affine types.  For each type, apply
`fixed_quadrilateral_unimodular_isolation.md` to its first four nodes.  Its
fifth node is at chord distance at most the length of the extracted window,
which by (5) is at most

```text
4C sqrt(R)/(N-4).                                     (8)
```

Let `c_*>0` be the minimum isolation constant over the finite set of types.
Under (1), `N-4` tends to infinity with `R`, so (8) is at most
`c_* sqrt(R)` for all sufficiently large `R`.  The fifth point then contradicts
the isolation theorem.  Hence (1) eventually fails, which is exactly (3).

Alternatively, the bounded-alphabet radius theorem in
`five_point_bounded_triangle_determinants.md` gives the same conclusion
directly, without passing through the finite list of affine types.

## 5. The precise obstruction to uniform extraction

For a window whose angular width is naturally `S`, equation (6) has scale

```text
R^2 S^3.
```

An endpoint arc has total angular width of order `R^(-1/2)`.  Averaging among
`N` points only guarantees a five-point window of width
`S=O_C(R^(-1/2)/N)`, and hence determinants of size

```text
O_C(sqrt(R)/N^3).                                     (9)
```

Therefore bounded unimodular type follows by this route only when
`N` is at least a constant multiple of `R^(1/6)`.  If `N` merely tends to
infinity along a sequence of radii, (9) may still tend to infinity.  The
fixed-quadrilateral theorem cannot be invoked because its constants depend on
the affine type, and (9) permits those types to escape to infinity.

This also explains why bounding endpoint triangle area by a multiple of
`sqrt(R)` is insufficient: determinant is invariant under `GL_2(Z)`, so a
quadrilateral of bounded coordinate height in a unimodular frame necessarily
has all triangle determinants bounded.  A new uniform argument must therefore
use arithmetic information beyond short-subarc averaging, or improve (9) by
an exact concyclicity mechanism coupling several separated windows.

Primitivity of the global difference lattice does not repair this loss.  For
example the four nodes

```text
(0,0), (q,0), (q+1,1), (0,1)
```

generate all of `Z^2` as a difference lattice: the triangle determinants
include `q` and `q+1`, whose gcd is one.  Yet one triangle determinant is `q`.
Since that determinant is unchanged by every unimodular coordinate change,
their coordinate height in every such frame is at least of order `sqrt(q)`.
This is a statement about general integer node sets, not an assertion that
these four nodes occur in an unbounded same-circle family.  It shows that
saturation or a primitive chord alone cannot bound moving affine shapes.
