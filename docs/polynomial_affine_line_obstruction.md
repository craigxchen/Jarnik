# No primitive affine-shape line for five polynomial circle points

This is a function-field theorem over `Q(i)[t]`. It excludes the
proposed degree-five polynomial circle family with quadratic
differences. It does not cover arbitrary sequences of integer
configurations, moving polynomial coefficients, or general rational
families without accounting for their denominators.

## 1. Statement

Let `m>=5` distinct polynomials `z_i(t) in Q(i)[t]` satisfy

\[
z_i\overline{z_i}=N(t),\qquad
\gcd_i z_i=1\quad\text{in }\mathbb Q(i)[t],
\]

where conjugation fixes `t`, and assume `N` is nonconstant. Define
the real polynomial triangle determinants

\[
D_{ijk}=\det(z_j-z_i,z_k-z_i),
\quad g=\gcd_{i<j<k}D_{ijk}\ \text{in }\mathbb Q[t],
\quad p_{ijk}=D_{ijk}/g.
\]

Then

\[
\boxed{\max_{i<j<k}\deg p_{ijk}\ge2.}
\tag{1}
\]

Thus the primitive Plucker map cannot have degree zero or one.
If `d=deg z_i` and `ell=max_(i<j) deg(z_i-z_j)`, then also

\[
\boxed{d\le3\ell-2.}
\tag{2}
\]

In particular `d=5` and `ell<=2` are impossible. Gaussian polynomial
primitivity can even be omitted from this particular corollary: a
common factor has degree at most two, and division by it leaves a
nonconstant primitive family to which (2) applies.

## 2. Primitive-triangle divisibility over the polynomial ring

Put

\[
b_m=\binom{\lfloor m/2\rfloor}3+
\binom{\lceil m/2\rceil}3.
\]

The local proof of the integer primitive-triangle theorem applies
over the two polynomial UFDs and gives

\[
\boxed{N^{b_m}\mid\prod_{i<j<k}p_{ijk}\quad\text{in }\mathbb Q[t].}
\tag{3}
\]

Here is the field-specific justification, including the possible
unsplit primes. An irreducible factor `P` of `N` over `Q[t]` either
remains irreducible over `Q(i)[t]`, or splits as `pi bar(pi)` up to
a constant. In the first case, conjugation preserves the prime, so
`v_P(N)=2v_P(z_i)` for every row. A positive valuation would give
a nonconstant common Gaussian polynomial factor, a contradiction.
Only the split case remains.

For a split factor with `P^e || N`, put `t_i=v_pi(z_i)`. Primitivity
at both conjugate primes gives `min t_i=0`, `max t_i=e`. The exact
triangle identity yields

\[
v_P(D_I)=e-(\max_{i\in I}t_i-\min_{i\in I}t_i)
+\sum_{\substack{\{a,b\}\subset I\\t_a=t_b}}r_{ab},
\]

where `r_ab=v_pi(z_a-z_b)-t_a>=0` for equal levels. This is exactly
the valuation formula in
[the primitive-triangle proof](affine_shape_radius_divisibility.md).
Its two cases, three or more valuation levels and exactly two levels,
give `sum_I(v_P(D_I)-min_J v_P(D_J))>=e b_m`. There is no exceptional
prime above two, since nonzero rational constants are units here.
This proves (3).

Every triangle polynomial is nonzero. Otherwise, for all sufficiently
large real `t`, three distinct points on a positive-radius circle
would be collinear. Distinct polynomial rows remain distinct outside
a finite set of parameter values.

## 3. Positivity excludes linear primitive triangle data

The polynomial `N` has no real root. If a real number `alpha` were a
root, writing `z_i=x_i+iy_i` would give

\[
x_i(\alpha)^2+y_i(\alpha)^2=0
\]

for every row. Thus every rational polynomial `x_i,y_i` vanishes at
`alpha`. Its minimal polynomial over `Q` would divide every `z_i`,
contrary to primitivity.

If every `p_I` had degree at most one, their nonzero product would
be constant or would have only real roots. Equation (3), with
`b_m>0`, would force the nonconstant polynomial `N` to divide that
product. This contradicts the preceding paragraph and proves (1).

For the degree consequence, the circumradius identity is the exact
polynomial identity

\[
4N D_{ijk}^2
=|z_i-z_j|^2|z_i-z_k|^2|z_j-z_k|^2.
\]

The squared norm of a nonzero Gaussian polynomial has twice its
degree: its leading real coefficient is the positive squared
modulus of its leading Gaussian coefficient. Consequently

\[
\deg D_{ijk}
=\deg(z_i-z_j)+\deg(z_i-z_k)+\deg(z_j-z_k)-d
\le3\ell-d.
\tag{4}
\]

Since some primitive determinant, and therefore some actual
determinant, has degree at least two, (2) follows.

## 4. Independent Gram-content proof for degree five

There is a second proof of the `d=5`, `ell<=2` exclusion. Some
difference must have degree two; otherwise (4) would give negative
triangle degree. Choose it as one side of a base triangle and let
`M` be the two-by-two Gram matrix of its two side polynomials.
Its maximal entry degree is four, while
`det M=D_base^2` has degree at most two by (4).

The primitive conic normal of the normalized five-point affine
shape has coefficient degrees at most three: its coefficients are
cubic products or sums of cubic products of the primitive triangle
polynomials, each of degree at most one. The exact formulas are in
[the feature-content audit](affine_conic_feature_content_audit.md).
Write `M=kappa Q`, where `Q` is the primitive polynomial Gram normal
and `kappa` is the full common polynomial content of the Gram entries.
The degree-four entry forces `deg kappa>=1`; the determinant identity
forces `deg kappa<=1`. Hence `kappa` is linear over `Q`.

At its real root, each diagonal Gram entry is a sum of two real
squares and vanishes. Both coordinates of both base sides must
therefore vanish. Thus `kappa` divides every side coordinate, making
`kappa^2` divide every Gram entry. This contradicts that `kappa` is
their full polynomial content.

## 5. Exact scope

The lower bound two in (1) is attained outside a short-arc setting.
For example, the five polynomials

\[
t+i,\quad -(t+i),\quad i(t+i),\quad t-i,\quad -(t-i)
\]

have common norm `t^2+1` and Gaussian polynomial gcd one. Their
triangle determinants have degree at most two, and include nonzero
constant multiples of both `t^2+1` and `t`; their common polynomial
gcd is therefore constant. This is not an endpoint construction.

The theorem excludes the proposed degree-five/quadratic-difference
family and every primitive degree-one Plucker family. It does not
exclude all higher-degree analogues: for instance `d=10`, `ell=4`
allows quadratic triangle polynomials, to which the real-root
obstruction above does not apply. No extraction of arbitrary integer
endpoint configurations into these polynomial hypotheses is asserted.
