# Moving elliptic maps: denominator growth with a logarithmic height margin

Fix an elliptic curve `E/Q`, a nontorsion point `P in E(Q)`, and a
positive presentation-degree bound `D`. There is an unconditional
moving-map denominator estimate when the joint coefficient heights
satisfy

```text
(H_n+1)(log(n+2))^6=o(n^2).                            (1)
```

For every nonconstant `f_n in Q(E)` with the bounded presentations
specified below, write its value in lowest terms as `a_n/b_n`, where
`b_n>0`. Then, for all sufficiently large `n`, this value is finite
and

```text
log b_n = (deg f_n) n^2 q(P)+o(n^2).                  (2)
```

Here `q` is the canonical height normalized by `h_x=2q+O_E(1)`.
In particular such nonconstant moving maps cannot have integral
values at `nP` for arbitrarily large `n`. Neither contact divisors nor
a determinant-one identity are assumed.

The logarithmic margin in (1) is part of the conclusion's hypotheses.
The argument below does **not** establish (2) under the weaker
condition `H_n=o(n^2)` alone. It also supplies no bounded-degree,
controlled-height presentation extracted from arbitrary integer
circle arcs, so it does not prove the uniform lattice-arc bound.

Section 8 records the additional fixed-splitting-field case. If all
denominator-section zeros lie in one fixed number field, Lang's fixed-target
theorem removes this logarithmic margin at the asymptotic `o(n^2)` level.

## 1. Presentation and the quantitative local assertion

Fix a nonsingular projective Weierstrass model. Write

```text
f=F/G,
F,G homogeneous of the same degree m<=D,
F,G in Q[X,Y,Z],    G|_E not identically zero,
H=h([coefficients F : coefficients G]).               (3)
```

The coefficients of numerator and denominator are charged jointly;
their relative scale is not discarded. Let

```text
L=log(n+2),             T=log(H+2),
R(n,H)=(H+1)(L+T)(1+log L+T)^5.                      (4)
```

There is a constant `C_(E,P,D)` such that, whenever `G(nP)!=0`,

```text
log max(1,|f(nP)|) <= C_(E,P,D) R(n,H).               (5)
```

The proof keeps the field degrees, target heights, coefficients of
the elliptic linear form, and local section constants uniform.

## 2. Primary elliptic-logarithm input

The primary source is Sinnou David,
[*Minorations de formes lineaires de logarithmes elliptiques*,
Memoires SMF 62 (1995), Theorem 2.1, printed p.10, and the first
remark on printed p.11](https://smf.emath.fr/system/files/filepdf/MSMF_1995_2_62__1_0.pdf).
Both pages were inspected directly, including the parameter formulas.

In David's notation, the nonzero linear-form bound contains

```text
(log B+log(d_K E_aux))
 (log log B+h_E+log(d_K E_aux))^(k+1)
 product_j log V_j,                                  (6)
```

times a constant determined by `k`, the field degree `d_K`, and
`log E_aux`. The initial condition `d_K log B>=log V_1` can be
removed: the p.11 remark adds `max_j max(0,log log V_j)` to each
of the two parenthesized factors in (6). This modification is
essential for moving targets; simply increasing `log B` to target
height would lose another factor of that height.

For a target `B_0 in E(Qbar)` of degree at most `M` and canonical
height `q(B_0)<=c(H+1)`, choose elliptic logarithms in a fixed
fundamental parallelogram. The distance from `nP` to `B_0` is
represented by

```text
Lambda=n u_P-u_(B_0)-r omega_1-s omega_2,
|r|+|s|<=C_E(n+1).                                   (7)
```

Apply the theorem with at most four copies of the fixed elliptic
curve. The points associated to these logarithms are `P,B_0,O,O`.
One may omit a zero logarithm. No independence of these points or
their logarithms is required for the theorem's nonzero-form bound.
For `nP!=B_0`, the representative (7) is nonzero.

The parameter choices and their dependence are:

* `K=Q(B_0)` contains the point coordinates and the integer linear
  coefficients; hence `d_K<=M`. Periods correspond to the rational
  origin and require no algebraic extension containing the periods.
* `h_E` is fixed. David uses the canonical height associated to the
  projective Weierstrass embedding; this is a fixed multiple of `q`.
* Choose `log V_(B_0)=C(H+1)` and all other `log V_j=C`, with a
  sufficiently large fixed `C`. Principal logarithms have bounded
  absolute value, so all the analytic requirements on `V_j` hold.
* Choose `E_aux=e`. Its upper-bound requirement holds after the
  same fixed enlargement of `C`.
* The integer coefficients in (7) have logarithmic height `O(L)`;
  take `log B=C L`. Use the p.11 modification instead of the
  discarded condition `d_K log B>=log V_1`.

Only one point-height factor grows. Thus
`product_j log V_j=O(H+1)` and
`max_j log log V_j=O(1+T)`. The exponent `k+1` is at most five.
The theorem consequently gives

```text
log |Lambda| >= -C_(E,P,M,c) R(n,H).                  (8)
```

With a fixed flat distance `dist_E` on the complex torus `E(C)`,
this gives

```text
log max(1,dist_E(nP,B_0)^(-1))
 <= C_(E,P,M,c) R(n,H),       nP!=B_0.                (9)
```

This derivation uses the 1995 theorem and its stated modification;
it does not presume an unverified improvement from a later theorem.

## 3. All denominator targets have bounded degree and height

The section `G|_E` of `O_E(m)` has a zero divisor of degree `3m`.
List its geometric zeros with multiplicity as `B_1,...,B_(3m)`.
Every individual target satisfies

```text
[Q(B_j):Q] <= 3D,
q(B_j) <= C_(E,D)(H+1).                              (10)
```

The first bound follows from the degree of the rational intersection
divisor. The second is the bounded-degree elimination/arithmetic
Bezout estimate proved in Section 2 of
[the uniform moving height lemma](moving_frame_uniform_height_lemma.md).
We apply (9) to one target at a time; no unbounded compositum of
fields is introduced. The constants depend only on `E,P,D`.

Under `H=o(n^2)`, none of these zeros equals `nP` for sufficiently
large `n`: otherwise (10) and `q(nP)=n^2q(P)` would contradict
nontorsion of `P`. This also handles zeros introduced by a
nonminimal presentation where both `F` and `G` vanish.

## 4. Uniform archimedean evaluation of the section ratio

Choose fixed smooth Hermitian metrics on the finitely many line
bundles `O_E(m)`, `1<=m<=D`, and fixed norms on their spaces of
sections. For a nonzero section `s` with zero divisor
`sum_j[B_j]`, one has the uniform estimate

```text
log ||s(X)|| >= log ||s||_section - C_(E,D)
                + sum_j log min(1,dist_E(X,B_j)).     (11)
```

Here section zeros are counted with multiplicity. To justify
uniformity, normalize the section norm to one. On fixed coordinate
discs, factor out the zeros using Weierstrass preparation. The
remaining nonvanishing factors vary continuously and stay bounded
away from zero on smaller discs, including as zeros coalesce when
the full multiplicity is factored out. The projective space of
normalized sections is compact. A finite cover of this parameter
space and of the curve gives (11). This is a statement about fixed
finite-dimensional spaces, not a moving-map constant left uncharged.

Scale the joint rational coefficient vector in (3) to primitive
integers. Its largest absolute coefficient is `exp H`. Therefore

```text
log ||F(X)|| <= H+C_(E,D).                            (12)
```

There is also a fixed lower bound on the section norm of nonzero
`G|_E`. Indeed restriction of integer polynomial coefficients to a
fixed rational basis of `H^0(E,O_E(m))` is a fixed rational linear
map. Its nonzero images lie in a fixed lattice, so their norms are
bounded below by a positive constant. Polynomial multiples of the
curve equation lie in the kernel and do not invalidate this claim.

Use this bound in (11) for `G`, and subtract from (12). Since
`F,G` are sections of the same metrized line bundle, their pointwise
metric factors cancel in the ratio. When `G(X)!=0`,

```text
log max(1,|f(X)|)
 <= H+C_(E,D)
       +sum_(j=1)^(3m) log max(1,dist_E(X,B_j)^(-1)). (13)
```

Equations (9), (10), and (13) prove (5). Degree-zero presentations
give constant functions and satisfy the local bound directly.

## 5. Denominator growth and the integral-value exclusion

Let `d_n=deg f_n`, the degree of its actual pole divisor after
cancellation, and set `Q_n=nP`. Apply the scalar case of
[the uniform moving height lemma](moving_frame_uniform_height_lemma.md),
equivalently to `[1:f_n:0:0:0]`, to obtain

```text
h(f_n(Q_n))
 = d_n n^2 q(P)
   +O_(E,P,D)(n sqrt(H_n+1)+H_n+1).                 (14)
```

Assume (1). Then eventually `H_n<=n^2`, so (4) gives

```text
R(n,H_n) <= C(H_n+1)(log(n+2))^6=o(n^2).             (15)
```

The global error in (14) is also `o(n^2)`. For a reduced rational
number `a/b`, `b>0`, there is the exact identity

```text
log b = h(a/b)-log max(1,|a/b|).                     (16)
```

Combining (5), (14)--(16) proves (2). A nonconstant function on an
elliptic curve has degree at least two. In particular these
denominators grow exponentially on the `n^2` logarithmic scale and
cannot equal one for all large `n`.

More generally, the same proof works whenever `H_n=o(n^2)` and
`R(n,H_n)=o(n^2)`; condition (1) is a simpler sufficient hypothesis.
The estimate uses the fixed archimedean place and an exact height
identity. It does not require a fixed set of finite primes.

For a fixed finite tuple `f_(1,n),...,f_(s,n)` of such moving functions,
the constants are uniform across its entries. Let `H_n` bound their
joint presentation heights, and suppose a positive integer `q_n`
clears every value. On the set of indices at which at least one entry
is nonconstant, its reduced denominator divides `q_n`, and hence

```text
log q_n >= 2n^2 q(P)+o(n^2).                         (16a)
```

The identity of that entry may vary with `n`: bounded degree, the
common height bound, and finite tuple size make the error uniform.
Thus subquadratic logarithmic clearing denominators are impossible
on any infinite set of nonconstant-frame indices. All-constant frames
require a separate argument using actual contact conditions.

## 6. What has and has not advanced

The new point beyond fixed-function quasi-integrality is the explicit
uniform dependence on the heights of moving poles and moving
numerator/denominator coefficients. The proof pays a logarithmic
loss and reaches the quantitative hypothesis (1).

No claim is made for all subquadratic heights: for instance
`H_n=n^2/log(n+2)` is subquadratic, but fails the sufficient condition
and is not controlled by (15). Removing this loss would require a
stronger argument. Applying the theorem to the lattice-circle
problem additionally requires an actual bounded-presentation-degree
extraction with the indicated coefficient-height margin; none is
provided here.

## 7. Arbitrary rational points on the fixed elliptic curve

The chosen orbit can be removed, at the cost of a rank-dependent
logarithmic exponent. Fix a basis `P_1,...,P_r` for the free quotient
of `E(Q)` and representatives for its finite rational torsion
subgroup. If `r=0`, there is no sequence of rational points with
unbounded height and the following statement is vacuous. Assume
`r>=1`.

Let `Q_n` be any sequence in `E(Q)` with

```text
h_n=q(Q_n) -> infinity,
```

and let `f_n` be nonconstant rational functions with presentations
(3) of degree at most the fixed `D` and joint coefficient heights
`H_n`. If

```text
(H_n+1)(log(h_n+2))^(r+5)=o(h_n),                    (17)
```

then `f_n(Q_n)` is finite for all sufficiently large `n` and its
positive reduced denominator `b_n` satisfies

```text
log b_n=(deg f_n)h_n+o(h_n).                         (18)
```

In particular integral values are excluded for all sufficiently
large indices. The exponent in (17) has not been replaced by the
rank-one exponent six when the rational Mordell--Weil rank is larger.

Here is the uniformity argument. Every `Q in E(Q)` can be written
uniquely in the form

```text
Q=T+sum_(j=1)^r m_j P_j,
T in E(Q)_tors,       m_j in Z.                      (19)
```

The fixed height-pairing matrix of this basis is positive definite.
If `lambda_min>0` denotes its smallest real eigenvalue, then

```text
sum_j m_j^2 <= q(Q)/lambda_min.                     (20)
```

The constants below may depend on the chosen basis, this positive
eigenvalue, and the finite torsion subgroup, all fixed once for
`E/Q`. There is no dependence on a newly chosen basis at index `n`.

For one denominator zero `B_0`, replace the moving target by
`B_0-T`. Rational torsion translation has a rational inverse, so

```text
Q(B_0-T)=Q(B_0),
q(B_0-T)=q(B_0)<=C_(E,D)(H+1).                       (21)
```

One may also take a maximum over the finitely many torsion
representatives for any ordinary-height or coordinate constants.
All these choices preserve the bounded target-field degree `3D`.
The elliptic distance to the target is represented by a nonzero
linear form, when `Q!=B_0`,

```text
Lambda=sum_(j=1)^r m_j u_(P_j)-u_(B_0-T)
       -a omega_1-b omega_2.                        (22)
```

Choose the target logarithm in the same fixed fundamental
parallelogram as before. The two integer period coefficients satisfy
`|a|+|b|=O_(E,basis)(1+sum_j|m_j|)`, hence are
`O_(E,basis)(sqrt(q(Q)+1))` by (20). All integer coefficients in (22)
therefore have logarithmic height `O_(E,basis)(log(q(Q)+2))`.

Apply David's theorem and the same p.11 modification with
`k=r+3`: the `r` fixed generator logarithms, one moving target,
and two fixed periods. The product of the point-height parameters
still has just one moving factor `O_(E,basis,D)(H+1)`. For
`h=q(Q)>=1`, put

```text
L_h=log(h+2),       T_H=log(H+2),
R_r(h,H)=(H+1)(L_h+T_H)(1+log L_h+T_H)^(r+4).
```

The resulting elliptic-logarithm bound, followed by the unchanged
local section estimate in Section 4, gives

```text
log max(1,|f(Q)|)<=C_(E,basis,D) R_r(h,H),
                                            G(Q)!=0. (23)
```

Under (17), one has `H_n=o(h_n)`. Thus the target height bound
(10) excludes `G_n(Q_n)=0` eventually. Also `H_n<=h_n` eventually,
and consequently

```text
R_r(h_n,H_n)
 <= C_r(H_n+1)(log(h_n+2))^(r+5)=o(h_n).             (24)
```

The uniform global height lemma now gives

```text
h(f_n(Q_n))
 = (deg f_n)h_n
   +O_(E,D)(sqrt(h_n(H_n+1))+H_n+1)
 = (deg f_n)h_n+o(h_n).                              (25)
```

Equations (16), (23)--(25) prove (18). No assumption that the
points lie in one cyclic subgroup was made. The stronger margin
(17), bounded presentation degree, fixed curve, and absence of an
arc-to-map extraction remain essential parts of the stated scope.

The common-clearing consequence also extends: for a fixed finite
tuple with common presentation bounds satisfying (17), on every
height-diverging subsequence with some nonconstant entry one has
`log q_n>=2h_n+o(h_n)`. Again the nonconstant entry may vary, and
all-constant frames are not excluded by this denominator argument.

## 8. Fixed splitting field: removing the logarithmic margin

We now return to the setting of Section 7. Let `K` be one fixed number
field, together with one fixed embedding `K -> C`. Suppose that every
geometric zero of every denominator section `G_n|_E`, counted with its
multiplicity, is a point `B_(j,n) in E(K)` under this embedding. There are
at most `3D` such zeros, since the presentation degree is at most `D`.
The field and embedding are fixed while the sections and their zeros may
vary with `n`.

Assume only

```text
H_n=o(h_n),                                           (26)
```

where `h_n=hat h(Q_n)` and `Q_n in E(Q)`. Bounded-degree arithmetic
Bezout gives, uniformly in `j,n`,

```text
hat h(B_(j,n)) <= C_(E,D)(H_n+1)=o(h_n).              (27)
```

The height-pairing inequality therefore gives

```text
hat h(Q_n-B_(j,n))
 = h_n+O(sqrt(h_n(H_n+1))+H_n+1)
 = h_n+o(h_n).                                       (28)
```

For each zero set `P_(j,n)=Q_n-B_(j,n) in E(K)`. Apply Lang's Theorem 2
at printed page 40 of [Integral points on curves](https://www.numdam.org/article/PMIHES_1960__6__27_0.pdf),
as detailed in [the fixed-function note](fixed_elliptic_frame_quasi_integrality.md),
over the fixed field `K` at the chosen archimedean place. Use the fixed
function `1/x` on `E`, whose zero at `O` has order two.
For every `rho>0`, the theorem says that the `K`-points satisfying

```text
|1/x(P)|_v <= c / H_E(P)^rho                       (29)
```

have bounded height. Since (28) tends to infinity, the points `P_(j,n)`
eventually lie outside this bounded exceptional set. Near `O`, the local
parameter distance satisfies `|1/x(P)|_v` comparable to
`dist_v(P,O)^2`; away from `O` the inverse distance is bounded. Hence, for
every `epsilon>0`,

```text
log+ dist_v(Q_n,B_(j,n))^(-1)
 <= epsilon h_n+O_(epsilon,E,K,D)(1)                 (30)
```

for all sufficiently large `n`, with the same threshold for every `j`.
Here choose rho sufficiently small relative to epsilon and the fixed
comparison constant `log H_E(P)<=C_E(hat h(P)+1)`, then use (28).
There are at most `3D` targets, so taking their maximum preserves (30).
If the finitely many archimedean embeddings of `K` are allowed to vary
among the evaluations, apply the same argument at each embedding and take
the maximum of the finitely many thresholds and constants.

Zeros with multiplicity cause no difficulty: the multiplicity is at most
`3D`, and multiplying (30) by that bounded multiplicity remains
`o(h_n)`. An eventual equality `Q_n=B_(j,n)` is impossible under (27) and
`h_n -> infinity`, because it would give `hat h(B_(j,n))=h_n`.

Insert (30) into the fixed Hermitian section estimate (13). The denominator
section has at most `3D` zeros, and therefore

```text
log+ |f_n(Q_n)| <= H_n+O_(E,D)(1)+3D o(h_n)=o(h_n).    (31)
```

The global height estimate (25), valid under `H_n=o(h_n)`, gives

```text
h(f_n(Q_n))=(deg f_n)h_n+o(h_n).                     (32)
```

Using the exact reduced-fraction identity (16), we obtain

```text
log b_n=(deg f_n)h_n+o(h_n).                          (33)
```

A nonconstant map from an elliptic curve to `P^1` has degree at least two.
Thus a common integer clearing the values of a fixed finite tuple of these
moving functions satisfies

```text
log q_n >= 2h_n+o(h_n)                                (34)
```

on every height-diverging subsequence on which at least one entry is
nonconstant. The identity of that entry may vary; bounded degree, bounded
tuple size, and the common fixed field make all errors uniform.

This is a corollary of Lang's classical fixed-target theorem after the
translation `Q_n -> Q_n-B_(j,n)`, not a new approximation theorem. It is
non-effective: the thresholds depend on `E`, `K`, the selected place, and
the chosen `epsilon`. If the splitting field varies with `n`, or if the
coefficient height is only `O(h_n)` rather than `o(h_n)`, this argument does
not remove the logarithmic margin. It also makes no arc-extraction or
global lattice-point claim.
