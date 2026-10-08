# Exact radius distortion for adaptive diagonal cotangent maps

This note isolates the one-parameter diagonal family behind the actual-point
double inversion. It gives an exact all-edge formula and a parameter-uniform
lower bound. The bound proves that averaging finitely many diagonal choices
cannot remove a close-pair raw-radius cost. It does not by itself obstruct
the normalized endpoint condition: when the parameter tends to infinity,
the same inequality permits the normalized target angle to decrease.

## 1. Exact pair formula

More generally, let `H_i=(x_i,y_i)` be distinct primitive integer
half-angle rows up to projective equivalence, and let M be a nonsingular
integer two-by-two matrix. Write `Y_i=M H_i`,
`D_(M,iell)=Y_i dot Y_ell`, and `E_(M,iell)=det(Y_i,Y_ell)`. Then

```text
E_(M,iell)=det(M) det(H_i,H_ell),
gamma_(M,iell)=gcd(D_(M,iell),E_(M,iell)),
n_(M,iell)=(D_(M,iell)^2+E_(M,iell)^2)/(gamma_(M,iell)^2 epsilon_(M,iell)),
```

where epsilon is 2 exactly when both reduced coordinates are odd. The exact
least squared radius of the image phases is the ordinary lcm of all
`n_(M,iell)`, by the [all-edge radius formula](integer_cotangent_lcm_height_target.md).
Thus this pair-vector formula handles an arbitrary rational
projective map after clearing its denominators; the diagonal calculations
below make the dependence on one adaptive parameter explicit.

Let

```text
H_0=(1,0),   H_i=(X_i,L),       1<=i<=k,
```

be half-angle rows for rational circle phases, with L>0 and distinct finite
integers X_i. For coprime positive integers a,d, apply the diagonal real map

```text
M_(a,d)=diag(a,d).
```

The projective map depends only on rho=a/d; coprimality fixes its primitive
integral representative. For a finite pair i<j, put

```text
D_ij(a,d)=a^2 X_i X_j+d^2 L^2,
E_ij(a,d)=a d L (X_j-X_i),
gamma_ij=gcd(D_ij,E_ij),
epsilon_ij=2 if D_ij/gamma_ij and E_ij/gamma_ij are both odd,
            1 otherwise.
```

The reduced cotangent vector of the relative image phase is

```text
(D_ij/gamma_ij, E_ij/gamma_ij).
```

Consequently its reduced Gaussian denominator norm is

```text
n_ij(a,d)=(D_ij(a,d)^2+E_ij(a,d)^2)/(gamma_ij^2 epsilon_ij).   (1)
```

For the anchor edge 0i, use D_0i=a^2 X_i, E_0i=a d L and the same formula.
The exact least squared target radius, including the transformed anchor, is

```text
N_(a,d)=lcm( n_0i(a,d), n_ij(a,d) : i<j ).                 (2)
```

No row primitiveness assumption is needed in (1): replacing a raw row
(aX_i,dL) by its coordinate-primitive row changes both entries of its pair
vector by the same factors, which are already included in gamma_ij. Formula
(2) is the ordinary all-edge lcm, so it retains opposite Gaussian
orientations and the prime above two.

The algebraic identity behind (1) is

```text
D_ij^2+E_ij^2
  =(a^2 X_i^2+d^2 L^2)(a^2 X_j^2+d^2 L^2).                  (3)
```

It is Lagrange's identity for the two transformed rows. The formula also
applies to rational diagonal parameters after clearing a common denominator.

## 2. A parameter-uniform lower bound

Because gamma_ij divides E_ij, and epsilon_ij<=2, (1) gives, for every
finite pair,

```text
n_ij(a,d)
 >= D_ij(a,d)^2/(2 E_ij(a,d)^2)
 = 1/2 * ((rho X_i X_j/L + L/rho)/(X_j-X_i))^2,
```

where rho=a/d>0. Therefore

```text
N_(a,d) >= 1/2 max_(i<j)
       ((rho X_i X_j/L + L/rho)/|X_j-X_i|)^2.              (4)
```

This is an exact arithmetic obstruction to hiding a diagonal distortion by
ordinary row gcds. It uses only the edge gcd, so additional prime or all-edge
lcm interaction can increase the radius.

For a fixed pair with `X_i X_j>0`, the real-valued expression on the
right is minimized at `rho=L/sqrt(X_i X_j)`. Its continuous minimum
therefore bounds every admissible rational parameter:

```text
inf_(rho>0 rational) N_(a,d)
 >= 2 X_i X_j/(X_j-X_i)^2.                               (5)
```

For any finite list rho_1,...,rho_h, applying AM-GM separately to each factor
in (4) gives the geometric-mean version

```text
( product_r N_(rho_r) )^(1/h)
 >= 2 X_i X_j/(X_j-X_i)^2                                (6)
```

for every positive pair X_i,X_j. Hence choosing or averaging a small set of
diagonal parameters cannot make all images cheap on a close positive
cotangent pair. The statement remains valid if each rho_r is required to be
rational and nonconformal (rho_r != 1); restricting the parameter set only
raises the infimum.

The anchor edge has the independent exact bound

```text
n_0i(a,d) >= (rho^2 X_i^2/L^2+1)/2.                       (7)
```

Indeed gamma_0i divides adL. This is useful when all finite pair differences
are large, while (4) is strongest on a close pair.

## 3. What the inequality does and does not imply at the endpoint

For positive X_i,X_j, the half-angle tangent of the transformed pair is

```text
|E_ij/D_ij| = |X_j-X_i| / (rho X_i X_j/L + L/rho).         (8)
```

Thus the same expression that forces the radius in (4) also measures angular
compression. If rho is very large, (4) grows quadratically in rho, while (8)
shrinks like 1/rho. After taking a fourth root, the lower bound supplied by
this pair alone decays like rho^(-1/2). This is why (4)--(6) are a genuine
radius-distortion lemma but do not close the uniform C sqrt(R) endpoint
problem: these lower bounds leave open whether angular compression can
compensate for the increase in the least radius.

For the double inversion in the companion note, a=3,d=1, so (1)--(4) apply
with rho=3. The large observed target lcm is therefore an all-edge
arithmetic effect; the present lemma supplies a closed-form lower certificate
for every pair, independent of numerical sampling.

There is a useful endpoint specialization when the chosen source anchor is an
endpoint. If all finite cotangents in that frame are positive and
`A=min_i X_i`, then the anchor edge at the minimizing row gives

```text
N_(a,d) >= rho^2 A^2/(2 L^2).                            (9)
```

The source phase width is `Theta=2 atan(L/A)`, while that same edge has
image width `2 atan(L/(rho A))`. Thus its normalized target contribution
is bounded below by

```text
2^(-1/4) sqrt(rho A/L) * 2 atan(L/(rho A)).              (10)
```

This lower bound tends to zero like `sqrt(L/(rho A))` as rho tends to
infinity. It quantifies why (9)--(10) do not settle uniformity; it does
not assert that the actual normalized target width tends to zero.

## 4. Source-anchor dependence and common scaling

There are two different operations that must not be conflated. Reanchoring the
already mapped target only replaces each target phase by a ratio with a target
phase; it leaves the full set of relative phases and the lcm (2) unchanged.
Reanchoring the source first and then applying the same displayed diagonal
matrix changes the matrix in physical coordinates, because a diagonal map does
not commute with a nontrivial rotation.

The exact source-anchor dependence is still elementary. For a source anchor
j, let `K_(i|j)=(x_(i|j),y_(i|j))` be a primitive integer half-angle
row representing `q_i/q_j`, with `K_(j|j)=(1,0)`. These rows are obtained
by primitive reduction of the Gaussian product `H_i conjugate(H_j)`.
Applying the diagonal map in this new frame gives, for every pair i,ell,

```text
D^(j)_iell=a^2 x_(i|j)x_(ell|j)+d^2 y_(i|j)y_(ell|j),
E^(j)_iell=a d (x_(i|j)y_(ell|j)-y_(i|j)x_(ell|j)),
```

and the exact squared radius `N^(j)_(a,d)` is the lcm of the corresponding
reduced norms from (1), including the anchor row `K_(j|j)`. Thus averaging
over source anchors is an average over genuinely different maps and can
change the radius. The general D,E lcm formula (1)--(2) applies in every
anchor frame. The specialized formulas (4)--(7) apply only after clearing
that frame's finite-row denominators to a common positive L; the positive-pair
versions (5)--(6) additionally require the two resulting finite cotangents to
have positive product. Without these hypotheses, use the general D,E formula
and do not invoke the positive-pair AM--GM conclusion.

Here is a small exact order-of-operations check. Start with rows

```text
H_0=(1,0), H_1=(2,1), H_2=(3,1),   a=2, d=1.
```

Applying the diagonal map first gives finite cotangents (4,6). The three
edge norms are `(17,37,629)`, so `N^(0)_(2,1)=629`. Reanchor first at
H_1: primitive reductions of `H_0 conjugate(H_1)` and
`H_2 conjugate(H_1)` give rows `(2,-1),(7,-1)`, which the same diagonal
map sends to `(4,-1),(14,-1)`, i.e. finite cotangents `(-4,-14)`. Their
edge norms are `(17,197,3349)`, so `N^(1)_(2,1)=3349`. This proves that
source-anchor choice is a real adaptive degree of freedom, while also showing
that the pair formula remains exact in every chosen frame.

A common rational multiplier of all target rows leaves every phase and every
reduced pair vector unchanged. Any successful adaptive argument must therefore
analyze the anchor-indexed family `N^(j)_(a,d)`, rather than invoke
anchor invariance before applying a nonconformal map.

The [exact checker](check_adaptive_diagonal_radius_distortion.py) passes
135 diagonal-map fixtures, 15 general matrix cases, and the source-anchor
example. It independently constructs the primitive Gaussian target by
Euclidean gcd/lcm arithmetic and checks all lower bounds with integers.
Root independently
checked the formulas and corrected the source-versus-target reanchoring
distinction before integration into the status note.
