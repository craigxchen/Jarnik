# Retaining the integer rank-one factors: exact cut content

This continues the [matching-residual audit](cut_matching_residual_rational_projection.md). Retaining
its rank-one factors gives exact selected-cut information, but the
factor-height improvement is confined to an explicit nested prime-power
overlap. In particular it gives no such improvement on squarefree
one-layer allocations. This does not exclude new inequalities coupling
these factors across many rows.

Assume a Gaussian-primitive tuple `z_i` has common odd norm `N`, and fix
an anchor `z_0`. A selected collection of gap layers gives
`z_i=beta_i w_i`, `Norm(beta_i)=Q`, and `Norm(w_i)=N/Q`.

## 1. The unreduced factor and its reconstruction congruences

Choose coprime integers `g_i,h_i` with

```
(z_i-z_0)/(z_i+z_0)=i g_i/h_i,
gamma_i=h_i-i g_i,
lambda_i=2z_0/gamma_i.
```

At the anchor take `(g,h)=(0,1)`; at an antipode take `(1,0)` and
interpret the ratio projectively. Then

```
z_i+z_0=h_i lambda_i,       z_i-z_0=i g_i lambda_i,
z_i=z_0 conjugate(gamma_i)/gamma_i.                 (1)
```

The Gaussian integer `lambda_i` is integral: an integer Bezout
combination of `h_i` and `g_i` expresses it as a combination of the
integral sum and the chord divided by i. Put

```
kappa_i=lambda_i conjugate(beta_i).
```

The exact rank-one factorization is

```
Qw_i-z_0 conjugate(beta_i)=i g_i kappa_i,
Qw_i+z_0 conjugate(beta_i)=h_i kappa_i.              (2)
```

Writing `kappa=k+i l`, the matrix in the previous note is
`[g,h]^T [-l,k]`. Its first integer factor is primitive. The content
of the second is `gcd(|k|,|l|)`.

The integral reconstruction conditions are exactly

```
gamma_i kappa_i = 2z_0 conjugate(beta_i),
conjugate(gamma_i) kappa_i = 2Q w_i,
Norm(kappa_i)=4QN/Norm(gamma_i).                    (3)
```

Conversely, the two Gaussian divisibilities by `2z_0` and `2Q` in (3),
together with the norm equality, reconstruct integral beta and w of
norms Q and N/Q and the original product in (1). The selected-cut
requirement further specifies their prime allocations, as follows.

## 2. Exact prime allocations and the common divisor

Fix `p^e || N`, choose a prime pi above p, and write

```
t_i=v_pi(z_i),       t_0=v_pi(z_0).
```

Primitivity gives rows at both levels zero and e. Let B(t) count the
selected gap layers at heights at most t. Then B is nondecreasing,
its integer increments are zero or one, `B(0)=0`, `B(e)=f=v_p(Q)`, and

```
v_pi(beta_i)=B(t_i),       v_bar(pi)(beta_i)=f-B(t_i).
```

Because gamma has coprime rational coordinates, pi and its conjugate
cannot both divide it. Taking valuations in (1) therefore gives

```
v_pi(lambda_i)=min(t_0,t_i),
v_bar(pi)(lambda_i)=e-max(t_0,t_i).
```

Consequently the retained factor has the exact allocations

```
v_pi(kappa_i)=min(t_0,t_i)+f-B(t_i),
v_bar(pi)(kappa_i)=e-max(t_0,t_i)+B(t_i).            (4)
```

There are no other odd prime factors: `lambda_i gamma_i=2z_0` and
beta has norm dividing N. Formula (4), together with unit choices,
is the complete source-prime information in the retained factor.

Let `d_Q=gcd_G(Q,z_0)`, up to a Gaussian unit. Then

```
gcd_G(kappa_i : all i) = (1+i)^r d_Q,              (5)
r=min_i v_(1+i)(lambda_i) in {1,2}.
```

For the pi valuation, the expression in (4) is nondecreasing for
`t<=t_0` and nonincreasing for `t>=t_0`. Its minimum is attained at
an endpoint and equals `min(f,t_0)`. The conjugate minimum is
`min(f,e-t_0)`. These are precisely the valuations of `d_Q`.
At the ramified prime, primitive gamma has valuation zero or one,
while `2z_0` has valuation two. Beta has odd norm. This proves (5),
including its entire two-part.

In particular the ordinary integer content of each kappa at an odd
split prime is the minimum of the two displayed exponents in (4).
There is no additional unexplained common divisor from the joint
rank-one representation.

## 3. The exact factor-height gain

Removing `d_Q` gives

```
Norm(kappa_i/d_Q)=Norm(lambda_i)/E_Q,
E_Q=Norm(d_Q)/Q
   =product_(p^e||N) p^min(f,e-f,t_0,e-t_0).        (6)
```

The identity follows from
`min(f,t_0)+min(f,e-t_0)-f=min(f,e-f,t_0,e-t_0)`.
The common ramified divisor from (5) is already present for lambda
and does not change this relative gain.

Thus all rows obtain the same height reduction, but only where both
selected and remaining layers occur at a prime and the anchor has
an intermediate allocation. If N is squarefree, `E_Q=1` for every
selection: the normalized kappa factors have **exactly** the original
lambda heights. Their phases can change from row to row, so their
original common short angular span does not automatically persist.
This calculation supplies no archimedean improvement in that regime.

## 4. Comparison with the primitive-chord divisor model

For a nonanchor row, let `c=gcd(|Re lambda|,|Im lambda|)` and
`s=sign(g)`. The primitive chord data of
[the independent attack note](endpoint_independent_attack.md), Section 1, are exactly

```
v=s i lambda/c,       d=|g|c,       t=-s h c,
z_i-z_0=d v,         2z_0=v(-d+i t).               (7)
```

Thus the primitive-chord divisor factorization is obtained by putting
the rational content c on the other factor. Equations (2)--(4) then
multiply its unnormalized lambda factor by the already chosen
Gaussian divisor `conjugate(beta_i)`. This is an explicit reversible
translation of the original divisor model plus the selected gap
allocations, not an extra per-row Diophantine equation.

For example, cross-row products satisfy exactly

```
kappa_i conjugate(kappa_j)
 =lambda_i conjugate(lambda_j) conjugate(beta_i) beta_j.  (8)
```

Their imaginary parts are integer determinants, but beta-dependent
rotations in (8) must be retained when giving an upper bound. Bounding
only the factor heights yields precisely (6); retaining the original
small lambda angle after these rotations would be an error.

A new joint constraint could still concern simultaneous residues or
cancellations among the factors in (8). This note establishes neither
their impossibility nor a sufficient bound. It identifies the exact
content and height bookkeeping that such a constraint must improve,
especially when the squarefree obstruction makes `E_Q=1`.

The [exact checker](check_cut_matching_rank_one_factor_equivalence.py)
verifies the allocation formulas, complete gcd, height identity,
reconstruction equations and primitive-chord dictionary on 2,304
anchored selections and 29,584 rows. These finite checks supplement
the general arguments above.
