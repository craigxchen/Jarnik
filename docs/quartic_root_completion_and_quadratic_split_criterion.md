# Exact root completion and quadratic-split criterion

This note records an exact reusable test for completing three prescribed
Gaussian-rational roots to a monic quartic with real cubic coefficients and
nonreal constant term, and an exact reduction for asking whether a real
rational quadratic can split all of its linear factors after composition.
It gives no global nonexistence result for the four-row problem.

## 1. Completing a quartic from three roots

Let \(r_0,u,v\in\mathbb Q(i)\), and put

\[
s_1=r_0+u+v,\qquad s_2=r_0u+r_0v+uv,\qquad s_3=r_0uv,
\]

with \(a=\operatorname{Re}s_1\), \(b=\operatorname{Im}s_1\),

\[
c=\operatorname{Re}s_2,\quad d_2=\operatorname{Im}s_2,\quad
e=\operatorname{Re}s_3,\quad f=\operatorname{Im}s_3.
\]

If the fourth root is \(r=x+iy\), reality of the first elementary
symmetric coefficient forces \(y=-b\). Reality of the second forces

\[
x=a-d_2/b.
\]

Thus, when \(b\ne0\), there is at most one completion, namely

\[
r=a-d_2/b-ib.
\]

The third elementary symmetric coefficient is real exactly when

\[
fb+a d_2b-d_2^2-b^2c=0. \tag{1}
\]

The constant coefficient is nonreal exactly when
\(\operatorname{Im}\(s_3r\)\ne0\). These conditions are necessary and
sufficient: the first three elementary symmetric coefficients are then
real, and the constant coefficient is \(s_3r\). All quantities are rational,
so a successful completion remains in \(\mathbb Q(i)\).

For a direct coordinate check of the completed set, write its four
roots as \(z_j=x_j+i y_j\), \(1\le j\le4\). The first three
elementary symmetric coefficients are real exactly when

\[
\sum_j y_j=0,\qquad
\sum_{j<k}(x_jy_k+x_ky_j)=0,\qquad
\sum_{j<k<\ell}(x_jx_ky_\ell+x_jx_\ell y_k+x_kx_\ell y_j-y_jy_ky_\ell)=0,
\]

with all sums over the four roots. For the completion formula above,
the first two equations hold by construction and the third is (1).
One separately checks \(\operatorname{Im}(e_4)\ne0\).

## 2. Exact test for a common real quadratic composition

Let the source roots of all seven factors be
\(z_1,\ldots,z_7\in\mathbb Q(i)\). A real rational quadratic
\(q(s)=A(s-h)^2+d\), with (A\in\mathbb Q^\times\) and (d,h\in\mathbb Q\),
has both roots of \(q(s)-z_j\) in \(\mathbb Q(i)\) if and only if

\[
z_j=d+A w_j^2\quad\text{for some }w_j\in\mathbb Q(i). \tag{2}
\]

The center \(h\) does not affect the square test: the two roots are
\(h\pm w_j\). Therefore all seven old factors split under one such
quadratic exactly when there are common \(d\in\mathbb Q\),
\(A\in\mathbb Q^\times\), and individual \(p_j,q_j\in\mathbb Q\) such that

\[
\operatorname{Re}z_j-d=A(p_j^2-q_j^2),\qquad
\operatorname{Im}z_j=2A p_jq_j\quad(1\le j\le7). \tag{3}
\]

This gives a reusable exact search condition. For a single root
\(z=x+iy\), (2) implies the necessary norm-square condition
\[
\sqrt{(x-d)^2+y^2}=|A|(p^2+q^2)\in\mathbb Q.
\]
It is not sufficient by itself; both equations in (3), or the exact
Gaussian-rational square test for \((z-d)/A\), must be checked. The common
\(d,A\) are essential: separate square tests with separate shifts or
scales do not establish a quadratic composition.

The condition is invariant under rational real affine changes of the
source variable \(z\mapsto \alpha z+\beta\), \(\alpha\ne0\): such a change
replaces \(d,A\) by \(\alpha d+\beta,\alpha A\). It is also invariant under
complex conjugation. These symmetries can be quotiented out in a search,
for example by fixing one source root's real part when it is distinct from
the others.

## 3. Audit against the known obstruction

For the seven roots in the balanced three-row family, two are
\(z_{13}=6-3i\) and \(z_{23}=6-4i\). If a common quadratic split them,
putting \(z=6-d\) in the norm-square consequences of (2) gives
\(z^2+9=u^2\) and \(z^2+16=v^2\) for rationals \(u,v\). The cited
quadratic-extension obstruction proves this forces \(z=0\), after which
the quotient of the two proposed Gaussian squares is \(3/4\), not a square
in \(\mathbb Q(i)\). Thus even this two-root necessary test rules out the
whole seven-root composition. The argument depends on these particular
imaginary magnitudes \(3,4\); it does not classify other triangles or
exclude other seven-root families.

Any finite search based on (3) can produce examples or suggest further
constraints, but failure in a bounded rational box is not a global
exclusion. A general classification would need to control the rational
parameters \(d,A,p_j,q_j\), as well as the row-completion equations.

## 4. A vertical-line triangle locus

There is a useful exact specialization of the three-row equations. Fix
the common root to i, and let the three pair roots be
X+ia, X+ib, X+ic, with X,a,b,c rational. For a row using the pair
X+ia, X+ib, the three prescribed roots are i, X+ia, X+ib. Writing
B=1+a+b, direct substitution gives

    s1=2X+iB
    s2=X^2-a-b-ab+iX(2+a+b)
    s3=-X(a+b)+i(X^2-ab)

Consequently its completion obstruction factors as

    C=(a+b)[X^2+(1+a+b)(1+a)(1+b)],

and, when B is nonzero, its private root is

    X(a+b)/(1+a+b) - i(1+a+b).

Assume X is nonzero, a,b,c are pairwise distinct, and a+b, a+c, b+c
are nonzero. Simultaneous completion of all three rows is then
equivalent to

    a+b+c=-2,    X^2=(1+a)(1+b)(1+c).

Indeed, the three nonzero pair sums allow us to set the bracket in C
to zero for each pair. None of 1+a,1+b,1+c can vanish because X is
nonzero. Subtracting the brackets for (a,b) and (a,c) gives
(b-c)(1+a)(2+a+b+c)=0. Thus a+b+c=-2, and any bracket now gives
the product identity. Conversely these two identities make every
bracket zero; their nonzero product also makes every completion
denominator nonzero.

Setting p=1+a, q=1+b, r=1+c gives the rational surface

    p+q+r=1,    pqr=X^2.

The seven source roots then have the compact form

    i;
    X+i(p-1), X+i(q-1), X+i(r-1);
    X(1+1/p)+ip, X(1+1/q)+iq, X(1+1/r)+ir.

The last three are the private roots opposite p, q, r, respectively.
The known triangle corresponds to (p,q,r,X)=(6,-2,-3,6). A sufficient
parametrization of points on the surface is p=s^2, q=-x^2, r=-y^2,
X=sxy, with s^2-x^2-y^2=1; this is a subfamily only, not a
classification.

The imaginary part of the completed row's constant coefficient, for the
row opposite r (whose pair roots have parameters p and q), is exactly

    X(r+1)(r-1)(r+pq)/r.

Thus the surface equation alone allows degenerate zero-constant cases;
for the desired nonzero constant one separately requires this expression
to be nonzero for each choice of missing parameter r, p, q.

For a modest exact search, I enumerated 44 distinct points on this surface
with rational p,q numerators of absolute value at most 8 and denominators
at most 4 (and r=1-p-q), then tested 420 common quadratic candidates per
surface point (18,480 total), derived from rational parameters t,y, each
with numerator and denominator bounded by 5 and 3, respectively. The
central root condition
i=d+K(ty+iy)^2 forces d=(1/t-t)/2 and K=1/(2ty^2), so each candidate
was checked by exact rational Gaussian-square tests against all seven
roots. The accompanying [checker](check_vertical_triangle_quadratic_split.py)
reproduces the finite enumeration. No all-seven split occurred in that
finite set. This is only search evidence: it neither excludes other
rational points on the surface nor other d,K.

## 5. A constructive rational family on the vertical-line surface

The surface is not just an isolated compatible triangle. For any rational
s,x,y with

    s^2 - x^2 - y^2 = 1,

set

    p=s^2,   q=-x^2,   r=-y^2,   X=sxy.

Then p+q+r=1 and pqr=X^2 identically. A two-parameter rational
subfamily is

    D=1-u^2-v^2,
    s=(1+u^2+v^2)/D,
    x=2u/D,
    y=2v/D,

for rational u,v with D nonzero. Substitution verifies the
hyperboloid identity exactly. This gives a genuine rational family of
constant-imaginary quartic triangles.

Within this chart, a simple exact nondegeneracy condition is

    D != 0, uv != 0, u^2 != v^2,
    2u != D, 2u != -D, 2v != D, 2v != -D.

The first two conditions make X and q,r nonzero; u^2 != v^2 makes
q != r. Also p=s^2 is positive and distinct from q,r. Since uv != 0,
the relation s^2-x^2-y^2=1 rules out p=1. The final four inequalities
are exactly x^2,y^2 != 1, so q,r != -1. Thus these conditions imply the
surface nondegeneracy conditions above. They exclude only finitely many
algebraic curves in the rational \(u,v\)-plane, so the admissible set is
infinite.

Here are exact nondegeneracy conditions. On the surface assume X nonzero,
p,q,r pairwise distinct, and p,q,r not in {0,1,-1}. The seven roots

    i;
    E_p=X+i(p-1), E_q=X+i(q-1), E_r=X+i(r-1);
    V_p=X(1+1/p)+ip, V_q=X(1+1/q)+iq, V_r=X(1+1/r)+ir

are distinct and nonreal, and no two are complex conjugates. Indeed, the
three E-roots have common real part and collide only when their
parameters agree; their conjugate pairs would require a pair summing to
2, hence the remaining parameter to be -1. The real parts X(1+1/t)
of the three V_t are injective in t. No V_t can equal or be conjugate to
an E-root because their real parts differ by X/t. The only possible
conjugate of the common root i among the private roots is V_{-1}=-i.
Nonreality of E_t,V_t excludes t=1,0, respectively. These conditions
also make all row constants nonzero.

For the row opposite parameter t, with the other two parameters a,b,
its four roots are i,E_a,E_b,V_t. Its polynomial is

    F_t(T) = T^4 - A1(t) T^3 + A2(t) T^2 - A3(t) T + A4(t) + i Y(t),

where

    A1(t) = X(3+1/t)
    A2(t) = 1+t^2+ab+3X^2
    A3(t) = X[1+t^2+(t+1)/t+ab(t^2-1)/t]
    A4(t) = X^2(3+1/t)+t^2
    Y(t) = X(t^2-1)(X^2+t^2)/t^2.

These formulas follow by multiplying the four linear factors and using
a+b+t=1 and abt=X^2. In particular all nonconstant coefficients are
real. Since

    Y(t)=X(t^2+X^2-1-X^2/t^2),

it is strictly increasing in t^2>0 when X>0, and strictly decreasing
when X<0. Therefore the three row imaginary constants are distinct
whenever p^2,q^2,r^2 are distinct. If two parameters were opposites,
their sum would be zero and the third would be 1, which is already
excluded. Thus the nondegeneracy assumptions imply distinct row
constants as well.

Rows opposite a,b share the full root i and the pair root E_c, where c
is the remaining parameter. For real T, their determinant polynomial
therefore has the exact factorization

    Im(conj(F_a(T)) F_b(T))
      = (Y(b)-Y(a)) (T^2+1) ((T-X)^2+(c-1)^2).

To see this, the left side is a real polynomial of degree four with
leading coefficient Y(b)-Y(a). It vanishes at i,E_c, as well as their
conjugates, because both row polynomials vanish at the two shared roots.
The displayed product has precisely those roots and the same leading
coefficient. Every determinant is consequently nonzero for real T.

For a small explicit member, take u=1/2,v=1/4. Then
s=21/11,x=16/11,y=8/11, and
(p,q,r,X)=(441/121,-256/121,-64/121,2688/1331). The checker
[check_vertical_triangle_family.py](check_vertical_triangle_family.py)
uses only Fraction arithmetic to verify the surface equations, all seven
distinct nonreal roots, absence of conjugate pairs, all three quartic
coefficient formulas, the three nonzero constant imaginary parts, and
the three determinant factorizations. This establishes a constructive
rational family; the finite quadratic-splitting scan above remains only
a separate bounded diagnostic.

Each nondegenerate rational member also supplies a full primitive
integer family, by the [polynomial specialization theorem](boolean_polynomial_integer_specialization.md).
This still has three nonanchor rows and four circle points; it does not
give a point count growing with the radius.
