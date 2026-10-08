# Full Boolean small-imaginary circuits and a split-norm gap

## Scope

The central full-Boolean core gives exact three- and four-row integer
circuits with **subpower coefficients** at the scale of each block.  Their
large factors have disjoint split-prime support.  This is a simultaneous
consequence of the small imaginary coordinates, stronger than applying
one core prime to one coefficient.  An explicit unbounded family shows
that one three-row circuit, even with positive split-supported factors,
bounded coefficients, the right logarithmic scale and a rank-two integer
determinant lift, need not yield a positive block-height gain.  The family
does not realize the prescribed Gaussian factors of the rows. A later
[unbounded three-row family](balanced_three_row_gaussian_polynomial_family.md)
does realize the full simultaneous system, with nonconstant fixed imaginary
coordinates. Thus a universal residual-height exclusion requires at least
four nonanchor rows. The larger systems and the radius-uniform arc count
remain open.

## 1. Exact circuits in the actual core

Fix `m>=3` nonanchor rows of the centrally truncated configuration and
write

```text
P_i=X_i+iY_i=K_i product_(nonempty T containing i) H_T,
n_T=Norm(H_T),
G_ij=product_(T containing i,j) n_T.
```

The blocks and their conjugates have disjoint rational prime support.
Assume their logarithms lie in `[(1-eta)w,(1+eta)w]`, where `eta=o(1)`,
and put `beta=max_i log max(1,|Y_i|)` and
`kappa=max_i log max(1,|K_i|)`.  The central reduction gives
`beta+kappa=o(w)` for fixed `m`.  The rows are distinct actual primitive
anchor numerators, so their pair determinants do not vanish.

Set

```text
D_ij=X_i Y_j-X_j Y_i=Im(bar(P_i) P_j),
t_ij=D_ij/G_ij in Z\{0}.
```

For each `T` containing `i,j`, the Gaussian integer `H_T` divides both
`P_i,P_j`.  Consequently `n_T` divides `Im(bar(P_i)P_j)`; disjoint prime
support makes their product `G_ij` divide it.  Moreover

```text
log |t_ij| <= log 2+beta+kappa+2^(m-1) eta w = o(w).       (1)
```

Indeed `|D_ij|<=2 max_k|P_k| max_k|Y_k|`, while
`log max_k|P_k|<=kappa+2^(m-2)(1+eta)w` and
`log G_ij>=2^(m-2)(1-eta)w`.  The bound is deliberately uniform over
all pairs.  It uses the small imaginary coordinates, not an equal-norm
claim for the `P_i`.

For three distinct rows put

```text
G_ijk=product_(T containing i,j,k) n_T,
E_jk;i=G_jk/G_ijk
      =product_(T containing j,k but not i) n_T.
```

The three `E` factors have pairwise disjoint rational prime support, and
each has logarithm `(2^(m-3)+O_m(eta))w`.  The rank-two minor identity
`Y_iD_jk-Y_jD_ik+Y_kD_ij=0` becomes the exact circuit

```text
Y_i E_jk;i t_jk-Y_j E_ik;j t_ik+Y_k E_ij;k t_ij=0.        (2)
```

Every multiplier in (2) has logarithmic size `o(w)`.  No division by a
Gaussian prime or appeal to a generic residue is involved.

When `m>=4` there is also a coefficient-only four-row form.  For
`I={i,j,k,l}` and a perfect matching `ij|kl`, define

```text
Q_ij|kl=product_(T: T intersect I is {i,j} or {k,l}) n_T.
```

The three `Q` factors are pairwise coprime, each with logarithm
`(2^(m-3)+O_m(eta))w`.  In every product `G_ab G_cd`, a block whose
intersection with `I` has size three contributes once, and one whose
intersection has size four contributes twice.  These contributions are
common to all three matchings; a size-two intersection contributes only
to its own matching.  Dividing the ordinary Pluecker identity by the
common product gives

```text
Q_ij|kl t_ij t_kl-Q_ik|jl t_ik t_jl
                   +Q_il|jk t_il t_jk=0.                 (3)
```

Its three coefficients also have logarithmic size `o(w)`.  Equations
(2)--(3) hold simultaneously, with the same `Y_i,t_ij,n_T`.
Related three- and four-row formulas appear for an older balanced
half-cut profile in
[the uniform-bound research note](codex_uniform_bound_research.md).
The point here is the full Boolean central core's subpower bound (1):
the large private factors in (2)--(3) each retain a fixed positive
multiple of `w` while their integer multipliers have height `o(w)`.

For the actual central-truncation output there is a sharper bound
that also tracks growing dimension. With `k=m+1`, original weight
`W`, removed half-mass `s`, and pair-moment error `zeta`, one has

```text
log|t_ij| <= log(C/2)+zeta W/4+2s.
```

The Gaussian gcd of the two primitive anchor rows has norm logarithm
`(d_0i+d_0j-d_ij)/2`; its excess over the retained common core is
at most the removed layer mass `2s`. Dividing the determinant by
this actual gcd norm cancels both anchor distances in its angular
bound. The [full-profile quantifier audit](endpoint_full_profile_quantifiers.md)
gives the proof and explicit constants. Individual residuals are
subpower when `(m+1)2^m/sqrt(M)->0`; products of all pair residuals
are subpower under the sufficient condition
`(m+1)^3 2^m/sqrt(M)->0`. Both include every fixed logarithmic
growth rate `m+1=floor(c log_2 M)` with `c<1/2`.

## 2. A sharp gap for one three-row circuit

For every even integer `a>4` with `5` not dividing `a`, let

```text
U=(a-3)^2+16,       V=(a+3)^2+16,       W=a^2+25.
```

Then

```text
U+V=2W,             log U,log V,log W=2 log a+o(1),      (4)
```

and `U,V,W` are pairwise coprime positive odd integers supported only on
split primes `p=1 mod 4`.  To see the support assertion, they are the
norms of the primitive Gaussian integers
`(a-3)+4i`, `(a+3)+4i`, and `a+5i`.  An inert prime dividing a
primitive Gaussian norm would divide both coordinates.  For coprimality,
any common odd prime of `U,W` or `V,W` divides `6a`; a prime dividing
`a` would then divide `25`, contrary to `5` not dividing `a`, while
`3` never divides `a^2+25`.  A common odd prime of `U,V` divides
`12a`, and the same argument applies.  Thus all three gcds equal one.

Assign the three private factors in (2) the values

```text
E_23;1=U,       E_13;2=W,       E_12;3=V,
Y_1=Y_3=1,     Y_2=2,          t_12=t_13=t_23=1.
```

Equation (2) is exactly `U-2W+V=0`.  The factors have full positive
height `2 log a+o(1)`, but every coefficient is bounded.

This circuit can even be lifted to honest rank-two integer minors at
the scale of a full three-row Boolean core.  Take
`a=2(5^h+1)` and `g=5^(2h)` for `h>=1`.  Then `g` is coprime to
`UVW`, is the norm of the conjugate-primitive Gaussian integer
`(2+i)^(2h)`, and `log g=2 log a+O(1)`.  Set

```text
(X_1,X_2,X_3)=g(3W,6W-V,2W),      (Y_1,Y_2,Y_3)=(1,2,1).
```

Direct calculation gives

```text
D_12=gV,             D_13=gW,             D_23=gU.       (5)
```

The common-triple block can also be made an **actual oriented Gaussian
divisor** without changing (5).  Write `(2+i)^(2h)=A+iB`, so
`g=A^2+B^2` and `gcd(B,g)=1`.  Choose the representative `r` of
`A B^(-1) mod g` with `|r|<=g/2`, and replace every `X_i` by
`X'_i=X_i+rY_i`.  The determinants are unchanged.  Since each old
`X_i` is a multiple of `g`, we have `X'_i=rY_i mod g`.  The identities

```text
Y_i A-X'_i B=0 mod g,
X'_i A+Y_i B=0 mod g
```

show directly that `(2+i)^(2h)` divides every `X'_i+iY_i` in the
Gaussian integers.  The conjugate does not divide them, because
`5` divides none of the `Y_i`.  Also `|rY_i|<=g`, much smaller than
the original `X_i`, so the logarithmic row scale is unchanged.

All `X_i` are positive, all three determinants in (5) are positive,
and therefore the slopes `Y_i/X_i` increase strictly with the row
index.  Also `log|X_i|=4 log a+O(1)`, while the three imaginary
coordinates remain bounded.  Thus (5) has the determinant
size, common-triple factor size and pair-private factor size predicted
by the formal three-row Boolean profile; the same conclusions hold for
the sheared `X'_i`.  The model still lacks the **pair-private** Gaussian
divisibilities of norms `U,V,W`.  At `h=1`, even after the common-block
shear `r=7`, the first row norm is `88 mod V=241`, so no Gaussian
block of norm `V` divides that row.  This missing simultaneous oriented
divisibility is exactly where a future full-system obstruction would
have to enter.

## 3. Completing the displayed pair orientations costs half a block

There is an exact height gap for completing this particular determinant
model with its **displayed** pair Gaussian factors.  Keep

```text
H_U=(a-3)+4i,       H_V=(a+3)+4i,       H_W=a+5i,
g=(a-2)^2/4,        (X_i,Y_i) as in (5).
```

Allow each entire `H_U,H_V,H_W` to be replaced by its conjugate.
Every integer rank-two lift with the same `Y=(1,2,1)` and the same
three determinants (5) has the form `X'_i=X_i+rY_i` for one integer
`r`: subtract two such lifts and use `Y_1=Y_3=1,Y_2=2` in their
three zero determinant differences.  Divisibility of the indicated
pair rows by the chosen orientations is equivalent to

```text
r = sigma_U (a-3)/4 - 2gW mod U,
r = sigma_V (a+3)/4 - 3gW mod V,
r = sigma_W a/5          mod W,                         (6)
```

where each `sigma` is `+1` or `-1`, and the fractions mean modular
inverses.  For example `H_U` divides `X'_3+iY_3` exactly when
`X'_3/Y_3` is `(a-3)/4 mod U`; the determinant `D_23=gU`
then transfers the same divisibility to row 2.  The other two
conditions work identically.  The denominators `4,5` are invertible
modulo their corresponding odd norms.

Let `M(a)=U(a)V(a)W(a)`.  Polynomial Chinese remaindering over
`Q[a]` gives a unique polynomial `R_sigma(a)` of degree at most five
whose residues modulo the three quadratic polynomials are the right
sides of (6).  Its degree-five coefficient is always positive:

| `sigma_U sigma_V sigma_W` | coefficient of `a^5` |
|:--|--:|
| `---` | `781/45000` |
| `--+` | `89/5000` |
| `-+-` | `311/18000` |
| `-++` | `319/18000` |
| `+--` | `311/18000` |
| `+-+` | `319/18000` |
| `++-` | `43/2500` |
| `+++` | `397/22500` |

The [checker](check_full_boolean_small_imaginary_circuits.py) computes
these remainders by exact rational polynomial arithmetic and verifies
all three defining congruences.  Every coefficient denominator of every
`R_sigma` divides `90000`.  Since `U,V,W` are pairwise coprime integers,
any integer solution of (6) satisfies

```text
M(a) divides 90000 (r-R_sigma(a)).                    (7)
```

Here the right side is an integer.  Indeed, after multiplying each
polynomial congruence by `90000`, both sides have integer coefficients;
division by each monic integer quadratic gives an integer quotient.
Now `R_sigma(a)=c_sigma a^5+O(a^4)` with `c_sigma>0`, whereas
`M(a)=a^6+O(a^5)`.  If `r=O(a^4)`, the integer in (7) is nonzero but
has absolute value less than `M(a)` for all sufficiently large `a`, a
contradiction.  More quantitatively, the same comparison gives
`|r| >= a^5/1000` for sufficiently large `a`, uniformly over the
eight signs.  Since the original `X_i` are `O(a^4)`, at least one
completed row has `log|X'_i+iY_i|>=5 log a+O(1)`.  The expected
four-block row height is `4 log a+O(1)`, so the displayed pair-block
completion costs at least `(1/2+o(1))w`, where `w=2 log a`.

This last theorem fixes the three primitive Gaussian factors shown
above, with only whole-block conjugation.  If a composite norm is
refactored and its rational primes are oriented independently, (6)
changes.  The theorem therefore does not rule out every Gaussian
factorization with the same three norms.  It does show exactly which
arithmetic information the bare rank-two circuit omitted.

One can vary the common norm in (5).  If it is a quadratic polynomial
`g(a)=g_2 a^2+g_1 a+g_0`, the same exact polynomial CRT calculation
gives degree-five coefficient

```text
c_sigma - g_0/300 + g_2/12,
c_sigma in {+/-13/90000, +/-3/10000, +/-1/4500}.       (8)
```

The coefficient is independent of `g_1`.  It cannot vanish when
`g(a)` has **integer polynomial coefficients**: multiplying (8) by
`300` would equate an integer with `300 c_sigma`, which is never an
integer.  For rational quadratic `g` the leading term can vanish
algebraically; for instance `g_0=0,g_2=9/2500` works with signs
`++-`.  Such a polynomial still has to take suitable positive
Gaussian-norm values, stay coprime to `U,V,W`, and admit a common
oriented block compatible with the CRT shear.  The calculation proves
none of those additional conditions, so it gives no full-system
counterexample or obstruction for moving `g`.

The countermodel does not settle the simultaneous four-row
relations (3), the near-square metric identities, or the full collection
of core-prime congruences.  It proves the narrower, unconditional point:
pairwise coprime split norm factors of the full expected height can obey
an isolated small-coefficient circuit, so the size and coprimality of
one such circuit alone give no fixed positive `w` gain.

The [exact checker](check_full_boolean_small_imaginary_circuits.py)
verifies twelve instances of (4)--(5), their coprimality and primitive
norm conditions, strict slope order, and the common-block Gaussian
divisibility.  The proof above supplies the unbounded family and the
polynomial CRT gap independently of those finite checks.
