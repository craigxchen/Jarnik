# Independent audit of the explicit balanced elliptic six-curve

This audit checks the construction in
[mixed_elliptic_six_curve_example.md](mixed_elliptic_six_curve_example.md)
independently of its intended application. The construction passes: it
gives a geometrically connected genus-one curve over `Q`, with a rational
origin, a birational map to `bar(M)_(0,6)`, degree one on every stable
boundary divisor, and positive rational rank. Exactly one forgetful image
is smooth; two more are singular rational anticanonical curves, and the
remaining three are birational images of anticanonical degree two.

This is a geometric example. It does not prove the global uniform
`C sqrt(R)` lattice-point arc bound or supply its missing arithmetic step.

## Cover and rational origin

For `A=1+2 lambda`, the equation of a coordinate is

```text
A x^2 -(A+lambda T)x-T=0,
Delta_lambda(T)=(A+lambda T)^2+4AT.
```

Its leading coefficient is the square `lambda^2`. The three monic
discriminants have root sets `{u,v}`, `{u,s}`, and `{v,t}` as claimed.
The squareclasses are independent even over `Qbar(T)`: the valuation at
`s` detects the second class, the valuation at `t` detects the third,
and either `u` or `v` then detects the first. Thus the normalization is
geometrically connected and has degree eight over the target line.

At each of the four distinct branch values the inertia has order two.
There is no branch at infinity, since each discriminant has even degree
and square leading coefficient. Riemann--Hurwitz gives
`2g-2=-16+4*4=0`. At `T=0` the independent roots are `0,1`, each simple,
so all eight points of this fiber are rational. Any one defines an
elliptic origin. The coordinate functions generate the entire extension,
so the normalized six-point map is birational onto its image. It extends
over the complete source by properness of the stable moduli space.

## Exhaustive boundary and contact audit

The homogeneous cross-product is exactly
`-(lambda-mu)X(X-Y)Y(X-2Y)`. It follows that equality of any two moving
coordinates occurs only at `0,1,infinity,2`. A collision with an anchor
also occurs only at the displayed special fibers.

At `T=0`, all eight binary choices occur. Each nonempty subset at zero
gives the cluster consisting of that subset and anchor zero, and each
nonempty subset at one gives the analogous cluster at one. This is
fourteen contacts, including both disjoint clusters at each mixed point.
At infinity, the three finite alternative poles are distinct and avoid
the anchors. Every nonempty subset choosing the infinite root therefore
gives precisely its anchored cluster, for seven contacts. At `T=2`,
the three alternative roots are distinct and avoid the anchors and `2`;
exactly the three moving pairs and the moving triple give contacts.

These twenty-five partitions are all distinct. More explicitly, their
partitions with a side of size two comprise the three anchor pairs,
nine anchor-moving pairs, and three moving pairs. Their partitions with
two sides of size three comprise the nine anchor-plus-two-moving
clusters and the moving triple. This is exactly `15+10=25` stable
boundary divisors.

All three target fibers are unramified. Independent differentiation gives
the claimed inverse derivatives

```text
0:        -1/7, -121/7, -25/7;
1:         4/7,   64/7,  16/7;
2:        7/15,  7/135, 7/39;
infinity: 7/3,  -7/57,  -7/9.
```

At an anchored cluster the anchor has slope zero, and all moving slopes
are distinct and nonzero. Blowing up that cluster once separates every
mark on its bubble. At the unanchored cluster at `2`, distinct moving
slopes similarly separate every mark after one blowup. Hence there are
no nested clusters, and each node has smoothing parameter a unit times
the source parameter. Simultaneous disjoint clusters at `T=0` each have
order one. Every listed pullback therefore has degree exactly one.

## Exact classification of all six projections

Forgetting `a`, `b`, or `c` is precisely the quotient by the corresponding
individual sign change. The `a` involution is absent from the inertia
list and is free; the `b` and `c` involutions each fix four points. Their
quotient normalizations thus have genus one, zero, and zero respectively.

For completeness, one can establish that the other three forgetful maps
are birational without the more involved general exclusion argument.
Boundary pullback gives `r delta=2` for every forgetful map, where `r`
is its degree onto the normalization and its image class is
`delta(-K)` on `bar(M)_(0,5)`. Thus any further nonbirational map would
have a distinct deck involution `sigma`.

For any pair of such deck involutions, the retained four-point cross-ratio
has degree four and is fixed by both. Distinctness follows because a
single involution associated with two labels would fix all normalized
coordinates. The pair therefore generates a Klein four group. For any
triple, using its complementary labels as the anchors shows that the
three involutions are independent: otherwise all three coordinates would
lie in one proper fixed field. Consequently `sigma` commutes with the
three known independent generators, and every triple among the four is
independent.

An elementary abelian two-subgroup of automorphisms of a characteristic-zero
elliptic curve has rank at most three: its translation subgroup lies in
`E[2]`, and its image in the origin-preserving automorphisms has rank at
most one. Thus `sigma=tau_a tau_b tau_c`. But `tau_a` is a translation
and `tau_b,tau_c` are reflections, so `sigma` is another translation.
The group generated by `tau_a,sigma` is a Klein four translation group.
Its quotient has genus one, whereas the degree-four cross-ratio forces
its fixed field to be rational. This contradiction excludes `sigma`.

The three known degree-two images have anticanonical class and arithmetic
genus one. The genus-one normalization forces the first image to be
smooth; the two rational normalizations force total singularity invariant
one for each of the other two. The remaining three images have class
`2(-K)`, arithmetic genus six, and genus-one normalization, hence total
singularity invariant five. They cannot be smooth.

## Quartic quotient, Weierstrass conversion, and positive rank

The kernel of the character `chi_b chi_c` is exactly the identity and the
three nontrivial free sign changes. It is the full translation group
`E[2]`, defined over `Q`. Multiplying the second and third discriminant
squareclasses gives the monic quartic

```text
F: Y^2=(T+7)(T+7/9)(T+7/81)(T+7/361).
```

There is no omitted quadratic twist: the discriminant leading factors
are rational squares. The quotient `E/E[2]` is isomorphic over `Q` to
`E` via the factorization of multiplication by two, with the quotient
origin understood. A different rational choice of origin on the quartic
is related by a rational translation and does not change the rank or
the torsion group up to isomorphism.

The explicit birational conversion to a cubic is

```text
x=3645-25200/(T+7),       y=230850 Y/(T+7)^2,
W: y^2=x(x+405)(x-35).
```

An exact verification, requiring no numerical approximation, is

```text
x     =3645(T+7/81)/(T+7),
x+405 =4050(T+7/9)/(T+7),
x-35  =3610(T+7/361)/(T+7),
3645*4050*3610=230850^2.
```

The rational point `(T,Y)=(0,49/513)` maps to `(45,450)`.
Both 13 and 17 are primes of good reduction: neither divides the
pairwise root differences `405,35,440` of the cubic. Exhaustive counting
of the cubic over the two fields independently gives
`#W(F_13)=20`, `#W(F_17)=24`. Prime-to-residue-characteristic torsion
injects on reduction. The second prime excludes any 13-primary torsion,
the first excludes any 17-primary torsion, and the remaining torsion
order divides both counts. The total rational torsion order therefore
divides four.

All four points of `W[2]` are rational, so the rational torsion is exactly
`(Z/2Z)^2`. In particular `(45,450)` is non-torsion because its
`y`-coordinate is nonzero. This proves positive rank for the source as
well. Independently, its twenty-four rational points over the three
unramified split fibers exceed the possible four torsion points.
An infinite non-torsion orbit meets the finite boundary support at most
finitely often, and hence does give infinitely many open configurations.

The supplied checker passed. This audit also separately evaluated the
cubic equation and exhaustively counted its finite-field points; the
boundary and projection arguments above check the geometric implications
that finite arithmetic assertions alone do not establish.
