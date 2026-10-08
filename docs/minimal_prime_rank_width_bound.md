# Bounded allocation widths at full prime rank

Fix `r>=1` and `C>0`. Let distinct Gaussian integers of modulus `R` lie
on an arc of length at most `C sqrt(R)`. Suppose there are exactly `r`
varying split-prime coordinates and their allocation rows have affine
rank `r`. Then every allocation width is bounded by an **effective
constant depending only on `r,C`**, uniformly over the primes, the
allocation matrix, the units, and the complete common Gaussian divisor.

If the sum of these widths is at least five, `R` is bounded by a constant
depending only on `r,C`. This second constant is generally ineffective.
In particular it applies at minimal support `r=M-1>=5` when the
allocation rows have full affine rank. That rank is automatic under
`C<=sqrt(2)` with arbitrary units, or `C<=2` in one literal unit class,
by [linear_allocation_affine_rigidity.md](linear_allocation_affine_rigidity.md).
It also applies when `r<=4` and the total width is at least five,
provided the stated rank hypothesis holds.

The new step is the effective width bound for a *moving integer allocation
matrix*. The fixed-matrix phase/Roth step is already in
[fixed_hadamard_roth_phase_obstruction.md](fixed_hadamard_roth_phase_obstruction.md).
The rank-one and weighted-power notes retain a phase kernel and hold the
integer matrix fixed. Here full affine rank removes that kernel first.
This hypothesis fails for the general extracted fair profiles, where
there are many more varying prime coordinates than independent rows.
No uniform bound on the number of points follows from this result.

## 1. Exact radius and phase coordinates

Remove the complete Gaussian gcd `d` of the rows, including all fixed
split allocations, ramified factors, and inert factors. Choose Gaussian
primes `pi_j` above the distinct varying rational primes `p_j`, with both
coordinates positive, and put

```text
z_i = d epsilon_i product_j pi_j^b_ij bar(pi_j)^(e_j-b_ij),
min_i b_ij=0,          max_i b_ij=e_j>=1,
N_0=product_j p_j^e_j,
D=log Norm(d)>=0,      W=log R^2=D+sum_j e_j log p_j,
H=max_j e_j,           E=sum_j e_j.                         (1)
```

The units `epsilon_i` are arbitrary. Select `r+1` affinely independent
rows, label one of them `0`, and form the nonsingular integer matrix

```text
B_ij=b_ij-b_0j          (1<=i,j<=r),
Delta=det B != 0,       J=adj(B),       L=r!.
```

The selected rows need not attain the extrema in (1). The elementary
Leibniz estimates still give

```text
|Delta| <= L H^r,
sum_i |J_ji| <= L H^(r-1).                              (2)
```

For `r=1`, use the usual empty determinant `1` in the adjugate.
Choose `phi_j=arg(pi_j)` in `(0,pi/2)` and lifts `y_i` of the point
arguments on the containing arc. With `x_i=y_i-y_0`, there are integer
vectors `t,l` such that

```text
2 B phi + (pi/2)t = 2 pi l + x,
|x_i| <= delta := C exp(-W/4).                           (3)
```

The common factor and the term `-sum_j e_j phi_j` cancel in the row
differences. Multiplying (3) by `2J` isolates every prime:

```text
4 Delta phi_j - pi m_j = 2 (Jx)_j,
m_j=4(Jl)_j-(Jt)_j in Z.                               (4)
```

Consequently, with `T=(2+2C)L`,

```text
|4 Delta phi_j-pi m_j|
    <= 2LC H^(r-1) exp(-W/4),
max(2|Delta|,|m_j|) <= T H^r.                           (5)
```

Indeed, `|m_j|<=2|Delta|+2LC H^(r-1)/pi` because `R>=1` and
`0<phi_j<pi/2`. Thus (5) has no dependence on the sizes of the unknown
integer phase lifts. In particular its coefficient bound does not
contain `W` or `log p_j`.

## 2. A uniform two-logarithm estimate

For an odd split prime `p`, set `alpha=pi_p/bar(pi_p)`, using
`log alpha=2i arg(pi_p)` and `log(-1)=i pi`. For integers `u,v` with
`u!=0`, put `B_*=max(|u|,|v|)`. Matveev's Corollary 2.3 gives the
convenient uniform specialization

```text
log |u log alpha-v log(-1)|
    > -K log p log(e B_*),             K=2^40.          (6)
```

Source: E. M. Matveev, *An explicit lower bound for a homogeneous rational
linear form in the logarithms of algebraic numbers. II*, Izvestiya:
Mathematics 64 (2000), 1217--1269,
[Corollary 2.3, printed page 1219](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&option_lang=eng&paperid=314&what=fullteng).
Its field degree is two; take `A_1=2 log p`, `A_2=4` and its
`C_1(2)<=2^32`. Its allowance to replace `B` by `B_*` gives (6),
because its numerical coefficient is at most
`2^32 * 2^2 * 2 * 4 * log(2e) = 2^32 * 32 * log(2e) < 2^40`.
Here `h(alpha)=(log p)/2` and both chosen logarithms have modulus at
most `pi`, so these height parameters are valid. The resulting numerical
constant is therefore below `2^40`.

For completeness, the minimal polynomial of `alpha` is
`p X^2-2(a^2-b^2)X+p`, where `pi_p=a+bi`; it is primitive and its roots
have modulus one, giving the stated absolute height. The linear form
in (6) is nonzero: otherwise exponentiation makes `alpha^u` a root of
unity, contradicting its nonzero valuation at `pi_p`.

Apply (6) with `u=2Delta`, `v=m_j`. Its linear form equals `i` times
the left side of (4). Combining (5) and (6) yields

```text
W/4 < K log p_j [1+log T+r log H]
          +(r-1)log H+log(2LC).                         (7)
```

Choose `j` with `e_j=H`. The **original** radius satisfies

```text
W=D+sum_h e_h log p_h >= D+H log p_j.                   (8)
```

In particular the complete common divisor strengthens (7); it has not
been discarded in the upper phase error. Since `p_j>=5`, (7)--(8) give

```text
H < A log(2H),
A=4 [ K (r+(1+log T)/log 2)
        +(r-1)/log 5
        +max(0,log(2LC))/(log 5 log 2) ].                (9)
```

This is uniform over both prime heights and all allocation matrices.
As `log(2H)<=2 sqrt(H)` for `H>=1`, an explicit coarse bound is

```text
H < 4 A^2,       H_0=ceil(4A^2).                        (10)
```

No limiting assumption that an allocation matrix stays fixed is used
in (7)--(10). The large constant is effective; it is not meant as a
practical search bound.

## 3. The finite-grid radius consequence

Once (10) holds, put

```text
D_0=L H_0^r,                 K_0=(L/2) H_0^(r-1),
T_0=union_(1<=h<=D_0) { pi m/(4h) mod 2pi : m in Z }.
```

This is one fixed finite set depending only on `r,C`. Equation (4)
implies, for every varying prime,

```text
dist(phi_j,T_0) <= K_0 C R^(-1/2).                      (11)
```

The fixed-target Gaussian phase lemma proved in Section 3 of
[fixed_hadamard_roth_phase_obstruction.md](fixed_hadamard_roth_phase_obstruction.md)
applies to each conjugate-coprime nonunit `pi_j`. It follows from
[Roth's theorem](https://www.cambridge.org/core/journals/mathematika/article/rational-approximations-to-algebraic-numbers/EFB89E2873019B246F004FAA06E05A7F)
in algebraic slope charts, with an elementary rational-slope case.
With epsilon `1/4`, it gives a fixed `c>0` such that

```text
dist(phi_j,T_0) >= c |pi_j|^(-9/4)=c p_j^(-9/8).
```

Thus, for a fixed `0<a<=1` depending only on `r,C`,

```text
p_j >= a R^(4/9).
```

These are bounds on the prime **norms**, not their Gaussian moduli.
Multiply them with the actual widths from (1):

```text
R^2=e^D product_j p_j^e_j
    >= e^D a^E R^(4E/9).
```

When `E>=5`, also `E<=rH_0`, so

```text
R^(2/9) <= e^(-D) a^(-rH_0),
R <= a^(-(9/2)rH_0).                                  (12)
```

This is a radius bound uniform over the allocation matrix and prime
support for each fixed `r,C`. It is generally ineffective because
the minimum of the finitely many Roth constants is ineffective.
When `E<=4`, this argument provides no strict radius-power gap.
The separate [five-point weighted-cube theorem](four_width_five_point_endpoint_gap.md)
handles that critical case directly: any five rows with `E<=4` have
`L>sqrt(2|d|)sqrt(R)`, or `L>2sqrt(|d|)sqrt(R)` in one common unit class.
Together the two results give bounded radius for every fixed
minimal-support `M>=5` in the small-arc ranges where affine rank is automatic.

## 4. Positive central weights at minimal support

Under `C<=sqrt(2)` with arbitrary units, or `C<=2` in one literal unit
class, [linear_allocation_affine_rigidity.md](linear_allocation_affine_rigidity.md)
gives strictly negative inner products for distinct rows

```text
f_i(j)=sqrt(e_j log p_j) (2b_ij/e_j-1).
```

Suppose now there are exactly `r+1` rows. They have affine rank `r`.
The centered integer rows `u_i(j)=2b_ij-e_j` consequently have rank
`r`, with one kernel relation on the rows. Its nonzero coefficients
cannot have both signs: the equality between the two positive sums
would have a strictly negative squared norm. They cannot include a
zero coefficient either, since pairing that omitted row with a positive
zero sum would give a strictly negative inner product. Reversing sign
and clearing denominators yields unique primitive integers

```text
n_i>0,       sum_i n_i u_i=0,       q=sum_i n_i,
sum_i (n_i/q)b_ij=e_j/2.
```

Cofactors give

```text
q <= (r+1) r! product_j e_j <= (r+1)r!H_0^r.            (13)
```

The exact product identity is

```text
product_i (z_i/d)^n_i = unit * N_0^(q/2).               (14)
```

Each `q e_j` is even because it equals `2 sum_i n_i b_ij`.
For lifts of the point arguments, their positive weighted mean lies
inside the same arc; relative to `arg(d)`, it belongs to
`(pi/(2q))Z`. Equations (10) and (13) therefore put this relative
central direction in a fixed finite torsion grid for fixed `r,C`.
The direction of `d` itself can vary with the radius; (14) does not
give a height bound on that absolute direction.

The pair-norm/obtuse package alone cannot bound `q`. For example, fix
the split primes `(1009,1013,1021)`, take any odd `e>=3`, and set

```text
b_0=((e-1)/2,(e-1)/2,(e-1)/2),
b_1=(0,e,e),       b_2=(e,0,e),       b_3=(e,e,0).
```

Their positive centered kernel is `(e,1,1,1)`, so `q=e+3` is unbounded.
Every pair satisfies the common-unit endpoint norm condition
`sum_j |b_ij-b_hj|log p_j > W/2+log 16` for `C=1/2` and `D=0`.
These allocations are realized by literal primitive Gaussian rows with
common unit, but their arguments are uncontrolled. The actual phase
condition, used in (3)--(7), is essential for bounding `q`.

## Verification scope

The [checker](check_minimal_prime_rank_width_bound.py) verifies exact
adjugate identities and coefficient bounds, actual Gaussian phase
inversion with arbitrary units and a common factor, the positive
cofactor relation, and the growing-denominator norm example. The
transcendence and Roth theorems are mathematical inputs, not verified
by this finite checker. No Lean formalization is claimed.
