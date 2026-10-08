# Complete Boolean fibers and the Gram-norm obstruction

The disjoint mixed elliptic six-curve has a useful exact divisor identity
that is slightly stronger than the count `7+7+7+4=25`.  It also makes
clear what elliptic multiplication can and cannot provide for the first
Gram lift.

## 1. The four complete fibers

Use the direction labels

```text
alpha       infinity       0       1       rho
T-value       16/3        -1/2     9/2      25/4
```

For an anchor `alpha` in `{infinity,0,1}`, let `P_(alpha,S)` be the
point in the fiber over its displayed `T`-value at which exactly the
moving labels in `S subset {a,b,c}` take the anchor direction.  Put

```text
p_alpha = P_(alpha,empty),
D_alpha = sum_(empty != S subset {a,b,c}) P_(alpha,S).
```

At `rho`, define `P_(rho,S)` in the same way, and put

```text
p_rho,empty = P_(rho,empty),
p_rho,j     = P_(rho,{j}),
D_rho       = sum_(|S| >= 2) P_(rho,S).
```

Every point here is rational, and every listed point has multiplicity
one in its unramified `T`-fiber.  If `F_tau` denotes the zero divisor of
`T-tau`, then the complete Boolean fiber decomposition is

```text
F_(16/3) = D_infinity + p_infinity,
F_(-1/2) = D_0         + p_0,
F_(9/2)  = D_1         + p_1,
F_(25/4) = D_rho       + p_rho,empty + p_rho,a + p_rho,b + p_rho,c.   (1)
```

The pole divisor of `T` has degree eight.  Therefore, for any two of
the four displayed values,

```text
div((T-tau)/(T-sigma)) = F_tau - F_sigma.              (2)
```

For example, the ratio between the `infinity` and `rho` fibers has
divisor

```text
D_infinity - D_rho
 + p_infinity - p_rho,empty - p_rho,a - p_rho,b - p_rho,c.
```

Thus the canonical fiber relation for grouped contact divisors retains
the private and empty cells.  The Boolean count by itself gives no
relation for the reduced `D`'s after those cells are deleted.

This identity has an elliptic multiplication interpretation.  Let
`pi:E -> E/E[2]` be the quotient by translations, identified with the
degree-four multiplication-by-two map after choosing the origin.  On
the quotient quartic, the two points above `T=tau` form a degree-two
divisor `Q_tau^+ + Q_tau^-`, and

```text
F_tau = pi^*(Q_tau^+ + Q_tau^-).
```

Pulling back the principal relation on the quotient gives (2).  The
multiplication map supplies complete fibers, not the reduced grouped
divisors with one or four cells deleted.

## 2. Coordinate sections recover the Boolean multiplicities

Let `X_j:E -> P^1` be the normalized coordinate for `j in {a,b,c}`.
For a direction `alpha`, let `H_(j,alpha)` be the degree-four zero
divisor of `X_j-alpha` (use `1/X_j` at infinity).  At an anchor fiber,
the direct Boolean multiplicity count gives

```text
sum_j H_(j,alpha)
 = D_alpha
 + sum_(|S|=2) P_(alpha,S)
 + 2 P_(alpha,{a,b,c}).                         (3)
```

At the nonanchor fiber it gives

```text
sum_j H_(j,rho)
 = D_rho
 + sum_j P_(rho,{j})
 + sum_(|S|=2) P_(rho,S)
 + 2 P_(rho,{a,b,c}).                            (4)
```

Equivalently, the right side of (4) is
`D_rho + (F_(25/4)-p_rho,empty) + P_(rho,{a,b,c})`.
Equations (3)--(4) are exact identities of effective divisors; they
explain why multiplying coordinate-difference sections produces
intersection multiplicities rather than the reduced contact section.

## 3. Why this does not produce the four Gram norms

The arithmetic divisibilities

```text
D_infinity | Norm(z),   D_0 | Norm(w),
D_1 | Norm(z+w),        D_rho | Norm(-3z+2w)
```

are congruences for an auxiliary integral Gaussian frame.  They are not
the assertion that the corresponding fixed rational functions on `E`
have zeros along the contact divisors.

There is a simple geometric obstruction to obtaining the required
simple contact zeros from fixed elliptic sections.  At every one of the
thirty-two Boolean fiber points, each moving coordinate is rational (the
inverse roots and the common targets are rational).  Hence a primitive
projective row `(p,q)` at such a point has

```text
p^2+q^2 != 0
```

over `Q`, including the directions `infinity,0,1,rho`.  A fixed Gaussian
linear section has norm divisor supported over the isotropic directions
`X_j=i` and `X_j=-i`, not over these real rational directions.

More generally, if `s` is a section over `Q(i)` and `P` is a rational
point of `E`, then

```text
ord_P(s * conjugate(s)) = 2 ord_P(s).
```

Thus a Gaussian norm section has even order at every rational point.
All contact points in (1) are simple.  In particular, a degree-eight
row-norm section cannot contain even one anchor divisor `D_alpha` of
degree seven: containing its seven rational points with norm order at
least two already costs degree `14>8`.  The analogous full fiber costs
degree at least `16`.  This rules out the desired low-degree fixed norm
sections.  It does not rule out a higher-degree norm section containing
`2D` (or additional zeros); such a section would not encode the
order-one contact divisibility used here.  The auxiliary frame norms
remain arithmetic values, and the determinant-one relation is a
pointwise Gram identity rather than a divisor identity on `E`.

The odd degree of `D_infinity,D_0,D_1` is a quick warning sign, but the
even-order argument is the real obstruction.  Adding `p_alpha` repairs
the degree count and gives the full fiber relation (1); it does not
repair the simple-point norm obstruction at degree eight.

## 4. The nonanchor private cells are optional

At `rho`, `H_(j,rho)` contains the singleton point `P_(rho,{j})`, but
the shared contact divisor `D_rho` contains only the three pair points
and the triple point.  A singleton point has only one moving coordinate
at `rho`, so no pair determinant is forced to vanish there.  Assigning
`P_(rho,{j})` to a private factor of row `j` would impose an additional
arithmetic congruence; it is not a consequence of the Boolean contact
profile.

The function `T-25/4` adds all four rho-private cells at once, while a
coordinate section singles out one moving label but also carries its
shared pair/triple zeros.  Neither construction forces the three
row-specific private factors.  This separates the necessary shared
divisibility from an optional private assignment.

## 5. Completed fiber forms and the exact residual equation

There is one useful arithmetic translation of (1).  If the rational
target is written in primitive coordinates `T=X/Y`, use the four linear
forms

```text
L_infinity = 3X - 16Y,
L_0         = 2X + Y,
L_1         = 2X - 9Y,
L_rho       = 4X - 25Y.                              (5)
```

They are fixed-denominator versions of `X/Y-tau`.  Choose the compact
recurrence subsequence in the integer-contact construction to avoid all
thirty-two Boolean points at the real place and at every prime in the
fixed bad-prime set.  Then the omitted local terms are bounded.  Passing
to a further subsequence fixes their valuations and the real signs, so
for fixed nonzero rational constants `c_alpha` the complete-fiber contact
contents satisfy `F_alpha = c_alpha |L_alpha(X,Y)|`.

Away from that fixed set, the valuation of `L_alpha` is exactly the sum
of the intersection depths with the eight points of the corresponding
complete fiber.  Consequently its contact content factors as

```text
F_infinity = D_infinity e_infinity,
F_0        = D_0         e_0,
F_1        = D_1         e_1,
F_rho      = D_rho       e_rho,                       (6)
```

where the constants `c_alpha` absorb the fixed-prime units.  Here
`e_infinity,e_0,e_1` are the missing anchor-cell contents, while `e_rho` contains the empty rho cell
and all three rho singleton contents.  This is the arithmetic form of
the complete fiber identities, and it does not choose a row for a rho
singleton.

Write the exact norm residuals as in the first-frame barrier and set

```text
x_infinity = R_infinity/e_infinity,
x_0         = R_0/e_0,
x_1         = R_1/e_1,
x_rho       = R_rho/e_rho.
```

Then the two Gram equations become the exact identities

```text
15 F_infinity x_infinity + 10 F_0 x_0
 - 6 F_1 x_1 - F_rho x_rho = 0,

4 F_infinity F_0 x_infinity x_0
 - (F_1 x_1 - F_infinity x_infinity - F_0 x_0)^2 = 4.   (7)
```

The private point height formula on this same recurrence subsequence
gives `log e_infinity,log e_0,log e_1=w_N+O(N)` and
`log e_rho=4w_N+O(N)`.  The size bounds in the first-frame barrier
therefore imply rational heights
`h(x_infinity),h(x_0),h(x_1) <= w_N+o(w_N)` and
`h(x_rho) <= 4w_N+o(w_N)` when `e` is retained in the denominator.
Since each `F_alpha` is a fixed linear form in `(X,Y)`, (7) is a genuine
low-height-coefficient quadratic/Pell-type system in the target ratio.
It is a useful exact accounting identity, but it does not force the
quadratic coefficients to vanish: the determinant-one equation is a
near-square equation with rational residuals, and rational-root bounds
apply only after a separate zero or fixed-coefficient argument.  Thus
elliptic multiplication supplies the four completed linear forms and
the accounting (6), while the existence or exclusion of the residual
Gram solution remains an arithmetic question.

The exact multiplicity check is in
[check_disjoint_elliptic_boolean_fiber_identities.py](check_disjoint_elliptic_boolean_fiber_identities.py).
