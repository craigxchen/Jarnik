# Exact conic compatibility and the quadratic-feature content

This audit uses the complete affine shape of a circle configuration.
It gives an exact six-point compatibility equation and computes the
primitive content and relation lattice of its quadratic-feature rows.
The compatibility is more information than the triangle-product
inequality alone. However, under the existing half-angle coordinates
it is an identity, and the feature content is an explicit function of
the same primitive triangle determinants. No new radius exponent gain
is proved here.

## 1. Normalized coordinates and their quadratic features

Let `m>=5` distinct integer circle points be labeled `0,...,m-1`.
Write `p_ijk` for the oriented triangle determinant divided by the
positive gcd of all triangle determinants; extend this notation
alternatingly to arbitrary ordered triples. Put

\[
D=p_{012},\quad U_i=p_{0i2},\quad V_i=p_{01i},
\quad W_i=D-U_i-V_i=p_{12i}.
\]

The affine coordinates relative to the first three points are
`(u_i,v_i)=(U_i/D,V_i/D)`. The first three become `(0,0),(1,0),(0,1)`.
For each remaining row define the integer quadratic feature

\[
F_i=\bigl(U_i(U_i-D),\ 2U_iV_i,\ V_i(V_i-D)\bigr).
\tag{1}
\]

If `d_1=z_1-z_0`, `d_2=z_2-z_0`, every row is orthogonal to

\[
\bigl(|d_1|^2,\ \operatorname{Re}(d_1\overline{d_2}),\ |d_2|^2\bigr).
\]

The feature matrix has rank exactly two. Indeed the circle becomes
a nonsingular conic through the normalized five-point subset. A
second independent conic through those five points would pull back
under its quadratic parametrization to a nonzero polynomial of
degree at most four with five distinct projective roots, which is
impossible. Thus those five points determine the conic uniquely.

Let `n=(a,b,c)` be its primitive integer normal, with sign chosen so
that `a,c>0`. Then

\[
a(u_i^2-u_i)+2b u_i v_i+c(v_i^2-v_i)=0,
\qquad ac-b^2>0.
\tag{2}
\]

Moreover `ac-b^2` is an integer square: the actual Gram determinant
is the square of the integer determinant of `d_1,d_2`, and division
by its common integer Gram content leaves a rational square that is
an integer. This condition retains the Euclidean realization; a
general rational conic need not have it.

## 2. All feature minors as cubic triangle products

For `3<=i<j<m`, put

\[
A_{ij}=V_iV_jp_{2ij},\qquad
C_{ij}=U_iU_jp_{1ij},\qquad
E_{ij}=W_iW_jp_{0ij}.
\]

Direct expansion gives the exact cross-product identity

\[
\boxed{F_i\mathbin\times F_j
=D\bigl(2A_{ij},\ A_{ij}+C_{ij}-E_{ij},\ 2C_{ij}\bigr).}
\tag{3}
\]

For example,

\[
Dp_{2ij}=U_iV_j-U_jV_i+D(U_j-U_i),
\]

so the first component of (3) is
`2V_iV_j[U_i(V_j-D)-U_j(V_i-D)]`, exactly the corresponding
feature minor. The third component uses
`Dp_{1ij}=U_iV_j-U_jV_i+D(V_i-V_j)`.
The first plus third minus twice the middle component equals
`2D W_iW_jp_{0ij}`, which verifies the remaining component.

Define `v_ij=(2A_ij,A_ij+C_ij-E_ij,2C_ij)`. Since all these vectors
are parallel to primitive `n`, there are uniquely determined integers
`h_ij` with `v_ij=h_ij n`. At least one is nonzero. Put

\[
h=\gcd_{i<j}|h_{ij}|
=\gcd_{i<j}\bigl(2A_{ij},A_{ij}+C_{ij}-E_{ij},2C_{ij}\bigr).
\]

The full gcd of all two-by-two minors of the integer feature matrix is

\[
\boxed{\operatorname{content}_2(F)=|D|h.}
\tag{4}
\]

Thus the common factor `D` in a naive quartic minor computation is
essential. After it is removed, every quantity in (4) is an explicit
cubic expression in primitive triangle data.

There is also an exact statement for the left relation lattice

\[
\Lambda_F=\{t\in\mathbb Z^{m-3}:\ \sum_{i=3}^{m-1}t_iF_i=0\}.
\]

Its rank is `m-5`, and its primitive exterior coordinates, up to signs
and complementation, are `h_ij/h`. Consequently

\[
\boxed{\det\Lambda_F=\left(\sum_{i<j}(h_{ij}/h)^2\right)^{1/2}.}
\tag{5}
\]

To see this without introducing a spurious metric factor, select any
two feature columns with a nonzero minor. Their minors are
`D h_ij` times one fixed nonzero component of `n`; dividing by their
own common gcd leaves precisely `h_ij/h`. The kernel is unchanged
by selecting these two independent columns. Dividing the sum of the
squares of *all three families* of feature minors by their joint gcd
would incorrectly insert an extra factor `||n||_2` into (5).

## 3. The exact six-point conic equation

For any six labels, written as `1,...,6` in this formula, the conic
condition gives

\[
\boxed{p_{123}p_{145}p_{246}p_{356}
=p_{124}p_{135}p_{236}p_{456}.}
\tag{6}
\]

This equality does not follow from the lower bound on the product of
all triangle determinants. But it becomes an identity under the
rational circle parametrization already used in the Gaussian
approach. To check this directly, fix a circle point `z_*` of norm
`N` and write every other point as

\[
z_i=z_*\frac{s_i+it_i}{s_i-it_i},\qquad
q_i=s_i^2+t_i^2,
\quad [ij]=s_it_j-s_jt_i,
\]

with rational projective parameters `(s_i:t_i)`. Then

\[
D_{ijk}=\frac{4N[ij][ik][jk]}{q_iq_jq_k}.
\tag{7}
\]

The two products in (6) contain exactly the same twelve pair brackets
after (7), and every row denominator occurs twice. Their constants
also agree. Dividing each triangle by its common integer content
preserves the equality. Thus neither side minus the other is a
nonzero arithmetic quantity to which a new size-versus-divisibility
argument can be applied. All higher conic feature determinants vanish
for the same reason.

This does not make (6) unnecessary for a proposed relaxation stated
only in triangle variables: such a relaxation must satisfy it. It
does mean that (6) adds no condition to a system already built from
actual half-angle rows by (7).

## 4. Radius dependence and the limitation

Let `H=max|p_ijk|`. Since every `U_i,V_i,W_i` is one of these
oriented triangle coordinates, (3) gives

\[
\|n\|_\infty\le3H^3,
\qquad
\det\Lambda_F\le
3\sqrt{\binom{m-3}{2}}\,H^3.
\tag{8}
\]

For the first bound, choose a nonzero `v_ij=h_ij n`; its entries
have absolute values at most `3H^3`, and `|h_ij|>=1`. For the second,
use a nonzero integer component of `n` to obtain
`|h_ij/h|<=3H^3` in (5).

On an arc of length `C sqrt(R)`, retaining the actual triangle gcd
`g` gives `H<=C^3 sqrt(R)/(8g)`. Hence (8) yields explicitly

\[
\det\Lambda_F\le
\frac{3}{512}\sqrt{\binom{m-3}{2}}
\frac{C^9R^{3/2}}{g^3}.
\tag{9}
\]

This is an upper bound, with no competing lower bound that supplies
a positive radius saving. In particular the cubic coefficient bound
for the conic does not establish `R=O(H^3)`: primitive Gaussian
factorization of its metric and the denominator of its center still
have to be controlled. The explicit cubic family shows that such a
universal exponent, if true, would be sharp, but does not prove it.

Equations (3)--(5) account for the complete feature-minor content and
its integer relation lattice. Equations (6)--(7) account for the
six-point conic identity. They neither rule out arithmetic uses of
several small actual values nor constrain a growing extracted cut
profile beyond what is proved here. The uniform endpoint bound
remains unproved.

The standard-library checker
[check_affine_conic_feature_content.py](check_affine_conic_feature_content.py)
independently verifies the cubic feature identities, primitive normal,
full minor content, and selected-column relation covolumes on 105 exact
integer circle configurations. It also checks every six-point identity
in those configurations. These checks supplement the displayed proofs.
