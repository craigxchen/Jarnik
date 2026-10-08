# A Kummer-cover bound for constant-imaginary Boolean polynomial rows

This note excludes **36 or more nonanchor rows** in a particular polynomial version of
the full Gaussian Boolean profile. It does not bound the number of lattice
points on an arbitrary circle: actual Gaussian blocks are integers, and an
integer profile need not interpolate to the constant-imaginary polynomial
family assumed here.

The geometric input is the Bogomolov–Miyaoka–Yau inequality for smooth
complex projective surfaces of general type. The section-arrangement setting
is described in [Urzúa, *Arrangements of rational sections over curves and the
varieties they define*](https://ems.press/journals/rlm/articles/5209),
Sections 2 and 7. We construct and calculate the double Kummer cover
directly. For the inequality, including its nonminimal form, see
[Hirzebruch, *Arrangements of Lines and Algebraic Surfaces*, §2.5](https://hirzebruch.mpim-bonn.mpg.de/244/1/hirzebrucharr.pdf)
or [Miyaoka's original theorem](https://doi.org/10.1007/BF01389789).

## Statement

Let `m>=3`, `d>=1`, and let `H_T in C[t]` for every nonempty
`T subset [m]` have degree `d`. Assume all `H_T` and all their
coefficientwise conjugates are pairwise coprime, including `H_T`
versus `bar(H_T)`. Roots within one block may repeat. Suppose nonzero
constants `K_i,c_i` and real
polynomials `X_i` satisfy

```text
P_i=K_i product_(T contains i) H_T=X_i+i c_i,
c_i in R\{0},                 X_i in R[t].                (1)
```

Then `m<=35` nonanchor rows. The conclusion is independent of `d` and of the
coefficients of the blocks. In particular, it applies to blocks in
`Q(i)[t]` and rows in `Q[t]+i Q`.

The coprimality hypothesis separates different block supports. No
assertion is made for varying imaginary polynomials `Y_i(t)` or
overlap between blocks or conjugates.

## Pair intersections and removal of the full common block

Set `f_i=X_i/c_i` and `D=2^(m-1)d`. Every `f_i` has degree at most `D`.
For each pair, the `2^(m-2)` blocks containing both labels give

```text
f_i-f_j=lambda_ij product_(T contains i,j) n_T,
n_T=H_T bar(H_T),         lambda_ij in R\{0}.          (2)
```

Indeed, the product on the right divides
`Im(bar(P_i)P_j)=c_i c_j(f_i-f_j)` and has degree `D`.
The difference is nonzero: otherwise the zero sets of `P_i` and `P_j`
would coincide, contrary to their distinct incident block sets.
Consequently equality holds up to a nonzero constant, and the leading
coefficients of the `f_i` are pairwise distinct.

The common full block gives `f_i congruent -i mod H_[m]` and
`f_i congruent +i mod bar(H_[m])`, including root multiplicities.
By the polynomial Chinese remainder theorem there is a unique
polynomial `h` of degree `<2d` satisfying both congruences.
Conjugating them shows `h in R[t]`. Put `F=n_[m]` and

```text
g_i=(f_i-h)/F.                                           (3)
```

Each `g_i` is real and polynomial. It has degree at most
`e=(2^(m-1)-2)d`, and its pair differences are

```text
g_i-g_j=lambda_ij product_(T contains i,j; T != [m]) n_T.
                                                               (4)
```

Thus each difference has degree exactly `e`; in particular the `g_i`
have distinct leading coefficients and no intersection at the base
point `t=infinity`. At each distinct root of `F`, (4) is nonzero,
so (3) introduces no new pair intersection.

For the next two sections assume additionally that each block is
squarefree. For every proper subset `T` of size `r>=2`, each of the
`2d` roots of `n_T` is exactly one
transverse `r`-fold intersection of the graphs of the `g_i` for `i in T`.
There are no other intersections. Therefore

```text
t_r=2d binom(m,r),                 2<=r<=m-1,          (5)
f_0=sum_(r>=2)t_r,
f_1=sum_(r>=2)r t_r,
t_2=d m(m-1).                                            (6)
```

The transformation (3) removes the full common block, even when its
roots repeat. Using it as a polynomial transformation avoids a
hypothesis that all original sections have empty total intersection.

## The smooth double Kummer cover in the squarefree case

Regard the graphs `S_i` of `g_i` as sections of the Hirzebruch
surface `F_e`, with negative section `C_0`, ruling fiber `F_0`, and

```text
C_0^2=-e,  C_0.F_0=1,  F_0^2=0,  S_i ~ C_0+e F_0.
```

Blow up once at every intersection of multiplicity `r>=3` and call
the resulting surface `Y`. Do not blow up the already normal-crossing
double points. Let `B` be the reduced union of the strict transforms
of the `m` sections and of those exceptional curves. It is a simple
normal-crossing divisor.

On the function field of `F_e`, with fiber coordinate `u`, adjoin

```text
w_i^2=(u-g_i)/(u-g_1),          i=2,...,m.              (7)
```

The divisor of each ratio is `S_i-S_1`; no ruling fiber or negative
section is a branch component. The valuations along `S_i`, `i>=2`,
prove that the `m-1` square classes are independent, so the cover
has degree `q=2^(m-1)`. Assign label `b_i` to the local double
monodromy at `S_i`; in `(F_2)^(m-1)`, the labels `b_2,...,b_m`
are a basis and `b_1=sum_(i=2)^m b_i`. Their only relation is
`sum_(i=1)^m b_i=0`.

At a proper `r`-fold point indexed by `T`, the exceptional curve
has label `b_T=sum_(i in T)b_i`, which is nonzero because
`empty != T != [m]`. At its intersection with a strict `S_i`,
`i in T`, the two labels `b_T,b_i` are distinct and hence
independent. The labels at an unblown double point are also
independent. Locally the normalization is therefore a product of
independent square-root coordinates and an étale factor. Its surface
`X` is smooth. Every component of `B` has ramification index two,
and at a node the two local inertia groups are independent. The
Hurwitz formula gives

```text
K_X=pi^*(K_Y+B/2).                                      (8)
```

In local coordinates this needs no covering formula beyond the
Jacobian: a smooth branch is `x=z^2`, so `dx wedge dy=2z dz wedge dy`;
at a node the independent labels give `x=z_1^2, y=z_2^2`, so
`dx wedge dy=4z_1z_2 dz_1 wedge dz_2`. Each branch therefore
contributes one half of its downstairs divisor to (8).

## Chern calculation

Write `J=K_Y+B/2`. On `F_e`, before the blow-ups,

```text
K_(F_e)+1/2 sum_i S_i
  ~ (m/2-2)C_0 + (e(m/2-1)-2)F_0.
```

At an `r`-fold blown-up point, the coefficient of the exceptional
curve in `J` relative to this pullback is `(3-r)/2`. Intersection
theory therefore yields

```text
4 J^2=m(m-4)e-8(m-4)
       -sum_(r=3)^(m-1)(r-3)^2 t_r.                     (9)
```

For the Euler number, stratify `Y` into the complement of `B`,
smooth parts of `B`, and nodes. The cover has respectively `q`,
`q/2`, and `q/4` points above a generic point of these strata.
There are `s=sum_(r>=3)t_r` blow-ups, so `e(Y)=4+s`;
`B` has `m+s` rational components, and it has
`t_2+sum_(r>=3)r t_r=f_1-t_2` nodes. Inclusion-exclusion gives

```text
4 c_2(X)/q =16-4m+f_1-t_2.                             (10)
```

By (8), `K_X^2=q J^2`. The pair-intersection identity

```text
sum_(r>=2) binom(r,2)t_r=binom(m,2)e
```

converts (9)--(10) into

```text
(3c_2(X)-K_X^2)/2^(m-3)
 =16-4m+3me+9f_0-2f_1-4t_2
 =16-4m+d[(18-m/2)2^m-4m^2-12m-36].                (11)
```

For every `m>=36` and `d>=1`, the last expression is strictly
negative: `18-m/2<=0`, and every other term is negative. The same
binomial sums in (9) give

```text
4K_X^2/q
 =d[(7m-36)2^(m-1)+m^2+3m+36]-8m+32>0
                                                     (m>=36). (12)
```

The positivity is immediate without a finite threshold check:
`(7m-36)2^(m-1)>0`, and
`d(m^2+3m+36)-8m+32 >= m^2-5m+68>0` for `m>=36`.

The cover is of general type in this range. Let `F` be a general
ruling fiber of `Y`. Its pullback is nef and, by (8),
`K_X . pi^*F=q(m-4)/2>0`. Hence for every integer `n>1`, the
divisor `(1-n)K_X` has negative intersection with a nef class and
has no nonzero section. Serre duality gives `h^2(nK_X)=0`.
Riemann--Roch and (12) now imply

```text
h^0(nK_X)>=chi(O_X)+n(n-1)K_X^2/2,
```

which grows quadratically. Thus `X` is a smooth projective surface
of general type. The Bogomolov--Miyaoka--Yau inequality
`K_X^2<=3c_2(X)` contradicts (11). It also applies if `X` is
nonminimal: each blow-up decreases `K^2` and increases `c_2` by one.
This proves `m<=35`.

The [dependency-free exact checker](check_boolean_section_kummer_uniform_bound.py)
passes 385 `(m,d)` substitutions (`m=4..80`, `d=1..5`) for the
squarefree intersection and Chern identities, and 2,018
proper-subset inertia checks through `m=10`. It also checks the
extensions below: 780 local repeated-root cases, 260 mixed
multiplicity profiles including `m=35,36,37,64`, and 84 base-genus
substitutions. These finite checks do not replace the all-`m`
calculation, the geometric resolution proof, or Bogomolov--Miyaoka--Yau.

## Repeated roots

The same bound holds without squarefreeness. Retain pairwise
coprimality of all blocks and conjugates, so a root belongs to one
`n_T` only. Equation (4) says that at a root of a proper `n_T`,
with `r=|T|>=2` and multiplicity `a>=1`, exactly the `r` sections
indexed by `T` meet. Every pair of them has contact order exactly
`a`; no other section passes through that point. Write
`b=a mod 2`, so `b` is zero or one. Sum below over all such
distinct roots, with their attached `r,a,b`.

Blow up the common point `a` successive times. After the `j`th
blow-up the new exceptional curve has inertia label
`(j mod 2)b_T`, where `b_T=sum_(i in T)b_i`; at each step all `r`
strict sections still have the same infinitely near point until the
last step, when their tangent directions separate. Thus precisely
the odd-index exceptional curves branch. Adjacent exceptional curves
have one branched and one unbranched label. If `a` is odd, the last
exceptional is branched and meets each strict section with an
independent label `b_T,b_i`; if `a` is even, it is unbranched and
the strict sections are mutually disjoint there. The normalized
Kummer cover of the resulting surface is smooth. We deliberately
blow up even an ordinary double point when `a=1`; this changes the
Chern bookkeeping but not the cover's validity.

For this local resolution, use the mutually orthogonal total
transform basis of the `a` exceptional curves. The coefficient of
the `j`th total transform in `K_Y+B/2`, relative to the pullback
from `F_e`, is `(3-r)/2` for odd `j` and `(1-r)/2` for even `j`.
The local subtraction from `4(K_Y+B/2)^2` is therefore

```text
ceil(a/2)(r-3)^2+floor(a/2)(r-1)^2.                 (13)
```

The `a` blow-ups increase `e(Y)` by `a`; the branch divisor gains
`ceil(a/2)` rational components; and only when `a` is odd do the
last exceptional and the `r` strict sections make `r` branch
nodes. The local contribution to `4c_2(X)/q`, beyond the base
constant `16-4m`, is exactly

```text
2a+(r-2)b.                                           (14)
```

Every proper block norm has degree `2d`, counted with multiplicity.
Consequently the weighted incidence totals are

```text
sum a = 2d(2^m-m-2),
sum r a = d m(2^m-4),
sum binom(r,2)a = binom(m,2)e,
e=(2^(m-1)-2)d.                                       (15)
```

Combining (13)--(15) with the unblown Hirzebruch intersection
number gives the **exact** repeated-root formulas

```text
4K_X^2/q
 =32-8m-3me+sum[(3r-5)a+(2r-4)b],

(3c_2(X)-K_X^2)/2^(m-3)
 =16-4m+3me+sum[(11-3r)a+(r-2)b].                    (16)
```

Because `0<=b<=a` and `r>=2`, (15)--(16) imply

```text
4K_X^2/q
 >=32-8m+d[(3m-20)2^(m-1)+4m+20]>0,

(3c_2(X)-K_X^2)/2^(m-3)
 <=16-4m+d[(36-m)2^(m-1)-16m-36]<0                 (17)
```

for every `m>=36,d>=1`. The first inequality uses
`2^(m-1)>=m`, giving the positive lower bound
`3m^2-24m+52`; the second is negative term by term.
The same general-type proof applies: the pullback of a general
ruling fiber is nef, `K_X` intersects it by `q(m-4)/2>0`, and
`K_X^2>0`. Riemann--Roch therefore gives quadratic pluricanonical
growth, while Bogomolov--Miyaoka--Yau contradicts the second line
of (17). This proves the statement without squarefreeness.

## Intrinsic base-genus corollary

Let `B` be a smooth projective complex curve of genus `g`, and let
`L` be a line bundle of degree `e=(2^(m-1)-2)d`, where `m>=36`
and `d>=1` are integers. Suppose `s_1,...,s_m` are distinct global
sections of `L`. For every proper subset `T subset [m]` with
`2<=|T|<=m-1`, let `D_T` be an effective divisor of degree `2d`.
Assume these divisors have pairwise disjoint supports and that the
zero divisor of every pair difference is exactly

```text
div_0(s_i-s_j)=sum_(proper T contains i,j) D_T.        (18)
```

Multiplicities within each `D_T` are allowed. Then

```text
d[(m-36)2^(m-1)+16m+36] <= 4(m-4)(g-1).              (19)
```

In particular, such arrangements with `m>=36` do not exist over a
curve of genus zero or one. The hypothesis concerns sections of one
line bundle, not arbitrary meromorphic functions with unaccounted
poles or intersections.

To prove this, compactify the total space of `L` to a ruled surface
over `B`. If `C_0` denotes its infinity section and `F_0` a fiber,
the graph sections have numerical classes

```text
C_0^2=-e,  S_i equivalent C_0+e F_0,
K equivalent -2C_0+(2g-2-e)F_0.
```

The functions `(u-s_i)/(u-s_1)` are globally defined ratios: the
local numerators transform by the same line-bundle transition
function. Their divisors are `S_i-S_1`, so the Kummer cover has
the same degree, branch labels, and smooth local resolutions as
above. Equation (18) ensures that all collisions are accounted for,
including at every point of the base.

The ruled surface has Euler number `4(1-g)`, and each original
section has Euler number `2(1-g)`. Exceptional curves remain
rational. Consequently the base constants `16-4m` and `32-8m`
in the repeated-root computation become respectively
`4(m-4)(g-1)` and `8(m-4)(g-1)`. With the same local notation
`r,a,b` as in (13)--(16), the exact formulas are

```text
4c_2(X)/q
 =4(m-4)(g-1)+sum[2a+(r-2)b],

4K_X^2/q
 =8(m-4)(g-1)-3me+sum[(3r-5)a+(2r-4)b].             (20)
```

The weighted incidence totals (15) are unchanged. Hence

```text
4K_X^2/q
 >=8(m-4)(g-1)+d[(3m-20)2^(m-1)+4m+20]>0,

(3c_2(X)-K_X^2)/2^(m-3)
 <=4(m-4)(g-1)-d[(m-36)2^(m-1)+16m+36].             (21)
```

The first lower bound is positive already at `g=0` by (17), and
increases with `g`. The pullback of a general ruling fiber is nef
and has intersection `q(m-4)/2>0` with `K_X`. Thus the same
Serre-duality and Riemann--Roch argument proves that `X` is of
general type. Bogomolov--Miyaoka--Yau makes the second quantity in
(21) nonnegative and proves (19).

Removal of a full common collision divisor is intrinsic as well.
If sections `f_i` of a line bundle `L_0` have a common collision
divisor `Z`, then `f_i-f_1` are sections of `L_0(-Z)`. Their pair
zero divisors are the original ones minus `Z`. Thus a full common
divisor of degree `2d` can be removed without a principal-divisor
assumption or a polynomial coordinate on the base.

This corollary is a function-field statement. A finite integer
circle configuration does not supply a family of sections over a
curve, and no arithmetic height bound follows from (19).

## Scope

The argument applies to a **constant-imaginary polynomial
interpolation** of every nonempty Boolean block, including repeated
roots within a block. It does not infer
such an interpolation from a large finite integer circle profile.
Nor does it rule out polynomial models with variable `Y_i(t)` or
overlap between different blocks or their conjugates. The proof gives no direct
radius-independent lattice-point bound.
