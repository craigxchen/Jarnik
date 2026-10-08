# Gauss reduction retains the marked nodes and their central weights

Reducing a positive binary Gram form is a valid symmetry of the six-node
description only if its marked projective vectors are changed at the same
time. The reduced form has bounded height at fixed determinant, but the
vectors can acquire unbounded height. The primitive central vector is
unchanged. This note records that distinction and corrects several tempting
small-determinant inferences; it proves no uniform endpoint bound.

## 1. An elementary reduction with its actual bound

Let

```text
Q = [[a,b],[b,c]],   gcd(a,b,c)=1,   Q positive definite,
det Q = m^2,        m a positive integer.
```

Choose a shortest nonzero vector of the integer lattice in the norm
`q(v)=v^T Q v`. It is primitive, since division by an integer greater than
one would shorten it. Extend it to an oriented integer basis, and shear
the second basis vector by an integer multiple of the first. The resulting
matrix `S` has determinant one and

```text
Q_red=S^T Q S=[[a',b'],[b',c']],
|2b'|<=a'<=c'.                                      (1)
```

The shear gives the first inequality and shortestness gives the second.
Consequently

```text
m^2=a'c'-b'^2 >= 3a'^2/4,
1<=a'<=2m/sqrt(3),
c'=m^2/a'+b'^2/a' <= m^2+m/(2sqrt(3)) <= 2m^2.       (2)
```

There are finitely many such integer triples at any fixed `m`. For `m=1`,
(2) forces `a'=1`, then (1) forces `b'=0`, and the determinant forces
`c'=1`. Thus a primitive positive determinant-one form reduces to the
identity. Entry-content and determinant are preserved by the unimodular
congruence.

These bounds concern `Q_red`, not the original canonical Gram matrix.
For example `[[1,n],[n,n^2+1]]` has determinant one and unbounded height.

## 2. Homogeneous central weights

Represent six distinct real projective nodes by nonzero columns
`v_i=(x_i,y_i)`. Write `d_ij=det(v_i,v_j)`. Up to a common scalar their
central weights are

```text
w_i = q(v_i)^2 / product_(j!=i) d_ij.                (3)
```

To verify (3), choose a real matrix `M` with `M^T M=Q`, write
`ell(v)=(Mv)_1+i(Mv)_2`, and set `z_i=conjugate(ell(v_i))/ell(v_i)`.
For `j=-2,-1,0,1,2`, the numerator of `w_i z_i^j` is

```text
ell(v_i)^(2-j) conjugate(ell(v_i))^(2+j),
```

a homogeneous quartic. The homogeneous Lagrange identity says that the
sum of its evaluations divided by the five determinants in (3) is zero.
Indeed, in a chart containing all six nodes this is the usual identity
`sum F(t_i)/product_(k!=i)(t_i-t_k)=0` for polynomials of degree at most
four. The five independent moment equations have a one-dimensional kernel,
as the same chart reduces their coefficient matrix to a Vandermonde matrix
times invertible factors. This proves the central-vector formula.

For rational nodes and rational `Q`, (3) is rational. Clearing denominators
and dividing by the ordinary gcd gives its primitive integer vector `u`,
well-defined up to sign; put `U=max_i |u_i|`.

Under a simultaneous change

```text
Q -> S^T Q S,       v_i -> S^(-1) v_i,       S in SL_2(Z),
```

every numerator and every determinant in (3) is unchanged. All raw weights
are therefore exactly unchanged. The matrix representing physical circle
points becomes `M S`, so `(M S)(S^(-1)v_i)=Mv_i`: these are the same points.

An independent change of representatives `v_i -> lambda_i v_i` multiplies
the numerator by `lambda_i^4` and its denominator by
`lambda_i^4 product_j lambda_j`. It consequently multiplies **every** raw
weight by the common factor

```text
(product_j lambda_j)^(-1).                           (4)
```

It does not introduce six independent rescalings. Primitive `u` and `U`
are invariant, including when a marked node moves to or from infinity.

In the chart with infinity first and five finite nodes `t`, (3) has the
opposite common sign to the convention

```text
w_infinity=-a^2,
w_t=(a t^2+2bt+c)^2 / product_(r!=t, r finite)(t-r).  (5)
```

This common sign is immaterial, but the infinity coordinate must be retained.

## 3. A complete small example

Take

```text
Q=[[2,1],[1,1]],
S=[[1,0],[-1,1]],
v=(infinity,0,1,2,3,5).
```

Then `S^T Q S=I`. Its inverse sends `(x,y)` to `(x,x+y)`, so the new
marked affine nodes, in the same label order, are

```text
(1,0,1/2,2/3,3/4,5/6).
```

Using (5), the original primitive central vector, with infinity first, is

```text
(-480,4,-375,3380,-6250,3721).                        (6)
```

The six new finite-node weights give exactly the same primitive vector up
to sign. Keeping the old affine nodes while replacing `Q` by the identity
would instead describe a different physical configuration.

The subgroup preserving all three distinct marked projective nodes
`infinity,0,1` is the identity in `PGL_2`. Integer shears fixing infinity
move the other two marks; restoring all three canonical marks restores the
original Gram point. Thus reduced height and canonical height cannot be
substituted for one another.

More quantitatively, if `v'=S^(-1)v`, (2) and positive definiteness imply

```text
q(v)=v'^T Q_red v' <= 4m^2 ||v'||_2^2,
||v'||_2 >= sqrt(q(v))/(2m).                         (7)
```

For canonical selected vectors infinity and zero, the two values on the
left are `a` and `c`. Large canonical entries at small `m` therefore force
large transformed integer vectors. Reduction moves the height into the
markings; it does not discard it. The constant in (7) is deliberately loose.

## 4. Explicit corrections to invalid small-determinant claims

The following statements are false and must not be used as lemmas.

* **A normalized node lies in `(1/L) Z` when `m=1`.** The actual clique
  `L=1`, finite numerators `1,2,3` has pair cotangents `3,2,7`, all
  integral. Selecting `A=1,B=3` gives `delta=2`, primitive Gram form
  `[[2,1],[1,1]]`, and the remaining normalized node `(2-1)/2=1/2`.
  The edge condition has the variable divisor `X-A`, not a fixed `delta`.
* **For a primitive form, `gcd(a,c)` divides `m^2`.** The primitive form
  `[[5,4],[4,5]]` has determinant nine but this gcd is five.
* **For a primitive form, `gcd(a,b)` divides `m`.** The primitive form
  `[[4,4],[4,5]]` has determinant four but this gcd is four.
* **An affine-chart change can itself lower `U`.** Equations (3)--(4)
  show the contrary when both the metric and every marking transform
  covariantly. A physical Mobius deformation with a fixed standard metric
  is a different operation, addressed in
  [the fixed-orbit audit](mobius_fixed_six_obstruction.md).

In addition, [the fixed-Gram content family](moving_base_fair_content_table.md)
has actual six-node cotangent cliques with a fixed determinant-one Gram
matrix and unbounded weighted evaluation content. It retains moving node
heights and does not contradict the fixed-base estimate. No general growth
improvement or uniform endpoint count follows from the reduction above.

The companion exact checker verifies the congruence, weight covariance,
independent rescaling, reduced bounds, and all displayed counterexamples.
