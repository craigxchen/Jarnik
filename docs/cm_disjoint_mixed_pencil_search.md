# Exact CM condition for the disjoint mixed quadratic pencil

This note records a concrete specialization test for the disjoint-support
balanced elliptic six-curve.  It is a bounded search result, not a global
nonexistence statement.

There is no obstruction at the level of the elliptic curve, rank, or
involution types.  For example

```text
E_5: y^2=x^3-25x
```

has `j=1728` and full rational two-torsion
`{O,(0,0),(5,0),(-5,0)}`.  The point `(25/4,75/8)` lies on `E_5` and is
non-torsion by the Lutz--Nagell theorem, since its coordinates are not
integral.  The three maps

```text
tau_a(Q)=Q+(0,0),   tau_b(Q)=-Q,   tau_c(Q)=(5,0)-Q
```

are independent commuting involutions: one translation and two
reflections.  They generate `E[2] semidirect {+1,-1}`.  The unresolved
condition is stronger: the three rational quotient coordinates must admit
four rational common arguments with distinct target values, so that they
belong to the required quadratic pencil and yield the twenty-five distinct
simple contacts.

Start with

```text
f_0(x)=x^2/(x-1)
```

and choose a monic split quartic

```text
R(x)=prod_j (x-r_j)=x^4+R_3x^3+R_2x^2+R_1x+R_0.
```

Put

```text
p_0=R_0,       p_1=R_1+R_0,
q_1=R_3,       q_0=R_2+p_1,
P_1=p_1x+p_0,  Q_1=x^2+q_1x+q_0.
```

Then, identically,

```text
x^2 Q_1-(p_1x+p_0)(x-1)=R(x).                 (1)
```

Thus all members

```text
f_lambda=(x^2+lambda P_1)/(x-1+lambda Q_1)     (2)
```

take the same value at the four rational arguments `r_j`.  Assume the
roots avoid all poles and branch fibers and give four distinct common
target values.

The member

```text
lambda_b=4p_0/p_1^2                             (3)
```

has branch value `0`.  Write

```text
u=p_1-4q_1,   v=p_0-4q_0,
lambda_c=(8u+4v-64)/(u^2+16v).                  (4)
```

The member `lambda_c` has branch value `4`.  For any parameter put

```text
D_1(lambda)=-2lambda p_1(1+lambda q_1)
            +4(-1+lambda q_0)+4lambda^2p_0,
D_2(lambda)=(1+lambda q_1)^2
            -4lambda(-1+lambda q_0).            (5)
```

The remaining branch values are exactly

```text
s=-D_1(lambda_b)/D_2(lambda_b),
t=-D_1(lambda_c)/D_2(lambda_c)-4.                (6)
```

Consequently the four branch points of the degree-eight Galois cover are
`0,4,s,t`.  Their cross-ratio, in the ordering `(0,4;s,t)`, is

```text
chi=(-s)(4-t)/((-t)(4-s)).                       (7)
```

The source has `j=1728` if and only if this unordered four-point set is
harmonic, equivalently if and only if one of the following three exact
equations holds:

```text
s t-2s-2t=0,
s t+4s-8t=0,
s t-8s+4t=0.                                    (8)
```

Indeed these are respectively `chi=-1,2,1/2`, the three harmonic
cross-ratios.  This gives an exact specialization equation, with no
dimension count.

The accompanying checker exhausts all four-element subsets of

```text
{-15,-14,...,15} minus {0,1,2}
```

and imposes the displayed necessary distinctness and denominator filters.
It finds no `j=1728` member even in that slightly oversized candidate set
of 20,452 quartics.  (The checker does not certify every member of this
superset as a noncancelled pencil with four simple contacts.)
This does not exclude rational common arguments.  After fixing three
arguments, substitution of (3)--(6) turns each equation in (8) into one
univariate rational equation for the fourth argument.  Factoring its
numerator and testing its rational linear factors is the remaining exact
finite search step for any chosen three-argument box.

The optional SymPy search
[search_cm_disjoint_mixed_pencil_rational_fourth.py](search_cm_disjoint_mixed_pencil_rational_fourth.py)
also fixes every three-element subset of that same integral box, factors
all three resulting univariate numerators, and tests every rational linear
factor.  It finds no rational fourth argument.  This remains a bounded
result because the first three arguments are restricted to the box.

Run
[check_cm_disjoint_mixed_pencil_search.py](check_cm_disjoint_mixed_pencil_search.py)
for the exact identity, harmonicity, and bounded-search checks.

## A rational harmonic component that fails the split-root condition

There is also a useful inverse parametrization of the first harmonic
case.  Write

```text
f_b(x)=c(x-l)^2/[c((x-l)^2-r(x-m)^2)].                 (9)
```

Assume `c!=0`, `r!=0`, and `l!=m`, so this is a degree-two map with
branch values `0,1`. Write the pencil in the affine-weight convention

```text
(P_k,Q_k)=(1-k)(P_0,Q_0)+k(P_b,Q_b)
```

and include its projective parameter at infinity. Eliminating the
parameter from the two branch-discriminant equations gives the
following necessary candidate condition:

```text
F(l,m,r)=
 27l^4-108l^3+108l^2
 +r(-60l^2m^2+168l^2m-96l^2
     +96lm^2-240lm+96l-24m^2+96m-96)
 +32r^2(m-2)^2(m^2-m+1)=0.                            (10)
```

This condition is independent of the representation scale `c`. It is
not sufficient at degenerate parameters: one must reconstruct a
rational common discriminant root, including a possible parameter at
infinity, and check that its numerator and denominator are coprime and define a
degree-two map. For example, at `c=1`, `(l,m,r)=(0,1,3/4)` satisfies
`F=0`, but its only projective common root gives a constant member.
On the nondegenerate locus these reconstruction checks recover the
desired member with branches `4,-2`.

The discriminant of `F` as a quadratic in `r` is

```text
144(m-2)^2 G(l,m),                                     (11)
```

where

```text
G=l^4m^2-16l^4m-8l^4+16l^3m^2+8l^3m+64l^3
  -12l^2m^2-24l^2m-48l^2-32lm^2+80lm-32l
  +4m^2-16m+16.
```

The component `l=2` is rational:

```text
r=9/[4(m^2-m+1)],       lambda_c=1/(3c+1).             (12)
```

The displayed affine parameter assumes `3c+1!=0`. For `c=-1/3` the
same pencil member lies at the projective parameter at infinity.

Thus harmonic branch values are not themselves an algebraic
obstruction.  On this component, however, the common-argument quartic
is

```text
x^2((x-2)^2-r(x-m)^2)-(x-2)^2(x-1),                   (13)
```

and its discriminant is

```text
-729(m-2)^4(53m^2-50m+50)
 / [2(m^2-m+1)^3].                                    (14)
```

For real `m != 2` this is negative.  Hence (13) cannot split into four
distinct rational roots, so this entire rational CM component cannot
produce the required twenty-five disjoint contacts.  Other rational
components of (10), and the other two harmonic orderings in (8), are
not excluded by this calculation.

## The rational twist on the first harmonic locus

The harmonic cross-ratio fixes `j`, but not the rational quadratic
twist.  In the parametrization (9), define

```text
N=-3(l-2)^2+4r(m-2)^2.                                (15)
```

Let `D_{2,b}` and `D_{2,c}` be the leading coefficients in `T` of the
two quadratic discriminants.  Direct calculation on the locus `F=0`
gives

```text
D_{2,b}=4c^2r(l-m)^2,

D_{2,c}=6c^2 r l^2(l-m)^2 N
          / [N+12cr(l-m)^2]^2.                        (16)
```

Therefore the product of their squareclasses is

```text
kappa_b kappa_c = 6N             in Q*/Q*2.            (17)
```

After deleting square factors, the elliptic quotient is a twist of

```text
Y^2=kappa T(T-4)(T-1)(T+2).                           (18)
```

The binary quartic on the right has invariants `I=108,J=0`.  In the
standard Jacobian convention, (18) has Jacobian

```text
y^2=x^3-2916 kappa^2 x.
```

Thus its congruent-number parameter `d=54kappa` has squareclass exactly
`N`.  The `j=1728` condition alone does not select a fixed rational
twist or guarantee positive rank.  On the rational component (12), the
twist parameter has squareclass

```text
N=4r(m-2)^2 ~ r ~ 1/(m^2-m+1).                        (19)
```

This twist calculation does not address the four-split-root condition;
equation (14) already excludes that condition on this component.

The [independent harmonic-component audit](cm_harmonic_component_audit.md)
records the exact resultant, its degenerate cases, and the parameter
at infinity. The polynomial equation is used only with those checks.

## The common-root quartic is not a natural two-isogeny quotient

One possible global obstruction would identify the elliptic curve
`y^2=R(x)` with a rational two-isogeny quotient of the source.  The known
disjoint example disproves such an identification.  Its source branch
cross-ratio and its common-root cross-ratio are

```text
lambda_E=-11090000009/27202483360,
lambda_R=5/3,
```

and `j(R)=438976/225`.  Starting from the Legendre model for `lambda_E`,
quotienting by each of its three rational two-torsion points gives three
`j`-invariants.  Modulo `101`, those invariants are respectively

```text
79, 1, 53,
```

whereas `j(R)` is `54`; all four denominators are invertible modulo
`101`.  Hence `y^2=R(x)` is not any of the three rational two-isogeny
quotients of the source.  The four rational diagonal points do not by
themselves furnish the proposed isogeny, so this route gives no global
CM obstruction.
