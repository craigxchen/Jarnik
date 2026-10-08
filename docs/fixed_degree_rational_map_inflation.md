# Rational maps of fixed degree and endpoint height

This note extends the small-projective-map calculation to rational maps
of any fixed degree. It uses one root in a finite extension of
`Q(i)` and keeps its ideal norms and all primitive-reduction costs.
It does not assume that the map splits over the Gaussian rationals.

The conclusion is a restriction on a possible descent, not the desired
uniform lattice-point bound. For `m` nonanchor rows, a nonexceptional
degree-`D` map of negligible coefficient height has

```text
log C_new >= (m-2D)W/8-o(W).
```

Thus, once `m>2D`, such a map cannot produce another endpoint cluster.
The exceptional pure powers are handled separately below. Maps of
growing degree or sufficiently large coefficient height remain open.

## 1. Input and normalization of the map

Let `h_i=x_i+i t_i`, `1<=i<=m`, be the primitive anchor numerators,
with

```text
x_i>0, t_i!=0, gcd_G(h_i,bar(h_i))=1,
|log|h_i|-W/4|<=E,
log|t_i|<=E, log|t_ij|<=E.
```

Here the exact primitive transitions are

```text
G_ij=gcd_G(h_i,h_j),
h_j bar(h_i)/N(G_ij)=x_ij+i t_ij,
x_i t_j-x_j t_i=N(G_ij)t_ij.                    (1)
```

Represent a nonconstant real rational map of degree `D>=1` by coprime
homogeneous integer forms `P(x,t),Q(x,t)` of degree `D`. Its circle
phase is `(P+iQ)/(P-iQ)`. Reanchor at the image of `(1,0)` by
multiplying `P+iQ` by the conjugate of its nonzero value there,
then clear any common scalar denominator. After this harmless
normalization write

```text
F=P+iQ,   F(1,0)=a>0,
H=max(1, absolute values of the Gaussian coefficients of F),
J=log H.
```

Thus `a` is a positive integer, `a<=H`, and `Q(1,0)=0`. The height
`J` in all statements is the height **after** reanchoring; this costs
at most twice the original logarithmic coefficient height plus a
degree-independent constant. The forms remain coprime.

Since `F` and `bar(F)` are coprime, a non-pure map has a root

```text
beta != i,-i,        F(beta,1)=0.                (2)
```

Indeed the coefficient of `x^D` is nonzero. If all roots were among
`i,-i`, having both would give a common factor with `bar(F)`.
The remaining cases are precisely `F=a(x+i t)^D` and
`F=a(x-i t)^D`.

## 2. One algebraic root gives almost coprime linear factors

Assume (2), put `alpha=a beta`, and work in
`K=Q(i,alpha)`, with `n=[K:Q]`. The usual monic rescaling of
`F(beta,1)=0` proves that `alpha` is an algebraic integer.
Cauchy's root bound gives, at every complex embedding of `K`,

```text
|sigma(alpha)| <= a+H <= 2H.
```

Define the algebraic integers

```text
L_i=a x_i-alpha t_i.
```

All gcds in this section mean ideals of `O_K`; divisibility retains
every prime-ideal valuation, including ramification. Let `g_ij`
be the ideal gcd of `(L_i)` and `(L_j)`. From (1),

```text
g_ij divides (a N(G_ij)t_ij).                   (3)
```

Modulo `G_ij`, the equality `x_i=-i t_i` gives
`L_i=-(alpha+i a)t_i`. The integer `t_i` is coprime to `G_ij`:
primitive integer coordinates give an integer Bezout identity for
`x_i,t_i`. Consequently

```text
gcd(g_ij,(G_ij)) divides (alpha+i a),
gcd(g_ij,(bar(G_ij))) divides (alpha-i a).
```

The two old Gaussian ideals are coprime, and remain so after
extension to `O_K`. Combining these facts with (3), valuation by
valuation, gives

```text
g_ij divides (a(alpha^2+a^2)t_ij).              (4)
```

The right-hand generator is nonzero by (2). Its absolute value
at each embedding is at most `5H^3 |t_ij|`. Therefore

```text
(1/n)log Norm(g_ij) <= E+3J+log 5.              (5)
```

Assume the explicit dominance condition

```text
W/4 >= 2E+J+log(8(D+1)).                       (6)
```

Then `x_i>=|h_i|/2`, and at every embedding the term `a x_i`
dominates `sigma(alpha)t_i`. In particular

```text
(1/n)log|Norm_(K/Q)(L_i)| >= W/4-E-log 4.       (7)
```

There is no unspecified field-degree factor in (5) or (7).

## 3. Exact denominator and primitive-content budget

Synthetic division over `O_K` gives

```text
L_i divides a^(D-1) F(x_i,t_i).                (8)
```

For example, if `F(x,t)=a x^2+b xt+c t^2`, then

```text
a F(x,t)=(a x-alpha t)(a x+(b+alpha)t).
```

The same division in degree `D` has integral coefficients because
`alpha` satisfies the monic rescaled equation.

Let `rho=Res(P,Q)` be the nonzero binary resultant. For coprime
integer `x_i,t_i`, the coordinate content

```text
c_i=gcd(P(x_i,t_i),Q(x_i,t_i))
```

divides `rho`. One proof uses the two homogeneous resultant Bezout
identities for `rho x^(2D-1)` and `rho t^(2D-1)`, and the coprimality
of `x_i,t_i`. The Sylvester determinant bound gives

```text
log|rho| <= 2D J+D log(D+1).                   (9)
```

After division by `c_i`, divide by `1+i` if both coordinates are odd.
The resulting Gaussian integer `d_i` is conjugate-primitive, and
the phase changes only by a Gaussian unit. Let `L` be a Gaussian
lcm of the `d_i`. Including the reanchored point with phase one,
the least possible circle radius is exactly

```text
R_new=|L|.
```

Indeed the integral anchor must be divisible by every `bar(d_i)`,
and `bar(L)` supplies an anchor attaining that radius. The units
introduced at two do not affect this denominator calculation.

Equations (8)--(9) show that **every** `L_i` divides `B L`, where
the one common scalar is

```text
B=a^(D-1)(1+i)rho,
log|B| <= (3D-1)J+D log(D+1)+(log 2)/2.         (10)
```

There is no product of row contents in (10). Since `B L` is
Gaussian, the norm of its extended ideal is `|B L|^n`.

For any list of integral ideals, sorting the valuations proves the
same lcm lower bound as over the Gaussian integers. Apply it to
the principal ideals `(L_i)`, using (5), (7), and (10). If
`K_m=binom(m,2)` and

```text
B_mD=m log 4+K_m log 5+D log(D+1)+(log 2)/2,
```

the result is the finite inequality

```text
log R_new >= mW/4-(m+K_m)E
             -(3K_m+3D-1)J-B_mD.               (11)
```

This argument does not require the polynomial to split over `Q(i)`
and does not assume pairwise coprimality without proving it.

## 4. A rational map cannot compress the angle beyond its degree

Put `y_i=t_i/x_i`. The polynomial `Q(1,y)` has a first nonzero
coefficient at some order `1<=e<=D`; that coefficient is a nonzero
integer. Under (6), its leading term dominates the others, while
`P(1,y_i)` remains positive. Thus

```text
|Q(1,y_i)/P(1,y_i)| >= |y_i|^e/[2(D+1)H],
|Q(1,y_i)/P(1,y_i)| <= 1.
```

The output angular separation from the anchor is twice its arctangent.
Since `|t_i|>=1`, `x_i<=exp(W/4+E)`, and `e<=D`, the angular span
of the output set satisfies

```text
log Delta_new >= -DW/4-DE-J-log(2(D+1)).         (12)
```

Combining (11)--(12) proves

```text
log C_new >= (m-2D)W/8
 -[D+(m+K_m)/2]E
 -[1+(3K_m+3D-1)/2]J
 -B_mD/2-log(2(D+1)),                           (13)
```

where `C_new=Delta_new sqrt(R_new)` is the normalized arc length.
In particular, for fixed `D` and `m>2D`, if `E+J=o(W)`, then
`C_new` tends to infinity. More generally (13) states every error
term when the row count also grows.

## 5. Pure powers and the scope of the result

The omitted maps are `F=a(x+i t)^D` and `a(x-i t)^D`. They raise
each relative phase to its `D`th power, with possible reflection.
If `R_intr` is the least radius of the input relative configuration,
their least output radius is exactly `R_intr^D`. On a sufficiently
small input arc their angular span is `D Delta`.

For a centrally extracted tuple with `m` nonanchor rows, the
independent nonempty core factors give

```text
log R_intr >= (1-2^(-m))W/2-o(W),
log Delta >= -W/4-o(W).
```

Consequently

```text
log C_new >= [D(1-2^(-m))-1]W/4-o(W).           (14)
```

For `D>=2` and `m>=2`, this again tends to infinity. The remaining
degree-one pure maps are the rotations and reflections already
identified in the projective note.

Thus no fixed-degree rational compression of negligible coefficient
height can supply a descent for arbitrarily large extracted clusters.
This does not address maps whose degree grows with the point count,
large coefficients outside (13), or transformations not given by one
common rational map. It proves no new bound on the original count.

## Verification

Two independent proof audits checked the algebraic ideal argument,
the binary resultant and common scalar, every coefficient in
(11)--(13), and the intrinsic-radius normalization in (14).
An exact Gaussian-arithmetic check also passed 200 degree-two and
degree-three triples, using seed 19407 and the maps

```text
2(x-(1+i)t)(x-(2+i)t),
2(x-(1+i)t)(x-(2+i)t)(x-(3+2i)t).
```

For primitive inputs with `3<=x<=100`, `1<=t<=20`, it checked
(4), (8), content divisibility by the binary resultant, the
primitive reduction, and `L_i | B lcm_G(d_i)`. These calculations
test the Gaussian-root specialization. Roots outside `Q(i)` are
covered by the audited ideal proof, not by this computation.
No Lean formalization is claimed.
