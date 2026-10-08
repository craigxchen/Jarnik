# A rank-four star pencil for the full four-row certificate

The [certified real algebraic four-row profile](four_row_real_polynomial_profile_certificate.md)
has a smaller exact representation: five monic quintics span a four-dimensional
real vector space, and all five arise from one complex polynomial pair of
degrees four and three applied to five real affine vector polynomials.
This representation is defined by rational operations on the certified roots.
It supplies necessary arithmetic conditions and a rational construction target;
it does not establish or refute rationality of those roots.

The dependency-free [checker](check_four_row_star_pencil.py) verifies the
nonzero minors and every division needed for this assertion by exact rational
interval arithmetic on the existing contraction cube. No decimal coefficient
is promoted to an exact algebraic relation.

## 1. The fifteen cuts as five vertices and ten edges

Use vertex labels `0,1,2,3,4`, with `0` the old anchor. Write `z_T` for the
fifteen roots in the original certificate, and set

```text
u_0 = conjugate(z_1234),            u_a = z_{a}       (a>0),
u_0a = conjugate(z_([4]\{a})),      u_ab = z_{ab}     (a,b>0).
```

These are the same fifteen cuts, consistently oriented toward the smaller
side of each partition of the five vertices. Define five monic quintics

```text
Q_a(t) = (t-u_a) product_(b != a) (t-u_ab),
n_ab(t) = (t-u_ab)(t-conjugate(u_ab)).
```

Put `Y_0=0` and retain the original four row constants `Y_a`. Directly
cancelling the common edge of `Q_0,Q_a` gives

```text
conjugate(Q_0) Q_a = n_0a P_a.
```

Thus, with `D_ab=Im(conjugate(Q_a)Q_b)`, one has `D_0a=Y_a n_0a`.
The elementary determinant identity

```text
Q_0 D_ab = Q_a D_0b - Q_b D_0a
```

then gives `deg D_ab<=2`: its right side has degree at most seven and
`Q_0` has degree five. Both `u_ab` and its conjugate are roots of `D_ab`.
The leading coefficient of the displayed identity is `Y_b-Y_a`, so

```text
D_ab = (Y_b-Y_a) n_ab.                                      (1)
```

All these constants are nonzero in the certified profile. In particular
all ten minors are quadratic and have negative discriminant.

This also makes the `S_5` symmetry transparent. After choosing a new
anchor `a`, all fifteen oriented roots undergo the same real affine map

```text
u -> (u-Re u_a)/(-Im u_a),
```

followed by a permutation of vertex labels. If the affine parameters are
`r=Re u_a`, `s=-Im u_a`, then the new quintics are
`s^(-5) Q_b(s t+r)`, with their labels permuted.

## 2. An elementary degree-two row reduction

Let `M` be the real two-by-five polynomial matrix with rows
`Re Q_a` and `Im Q_a`. Its minors are (1). Suppose its two row degrees
are `r,s` and `r+s>2`. Their leading coefficient vectors must be
proportional, since their exterior product is the coefficient of
`t^(r+s)` in the minors and all these coefficients vanish.
Subtracting a suitable real scalar multiple of `t^|r-s|` times the
lower-degree row from the other row therefore lowers a row degree.
This is an elementary polynomial row operation of determinant one.

Repeating terminates with row degrees summing to two. The sum cannot be
smaller, since at least one minor has degree exactly two. Thus there is
`U in SL_2(K[t])`, where `K` is the real coefficient field, such that

```text
U M = [ A_a(t) ]
      [ B_a(t) ],       deg(A)+deg(B)=2.                    (2)
```

No coprimality assumption on the ten minors is needed for this argument.
The possibilities are `(0,2)`, `(1,1)`, or `(2,0)`.

There are at most four coefficient vectors among the two rows in (2).
Since polynomial row operations do not change constant linear relations
between columns, the five `Q_a` satisfy a nonzero relation

```text
sum_a lambda_a Q_a = 0,       lambda_a in K.                (3)
```

Their real constant span consequently has dimension at most four.
This is a universal consequence of the full degree-one profile, not an
observed approximate relation peculiar to the stored root.

## 3. Exact additional information for the certified root

For the certificate, take the four real coefficient rows of the quintics
in degrees `5,4,3,2`, forming a four-by-five matrix `C`. The signed
maximal minors of `C` are approximately

```text
(-0.00284059700440246,
  0.01075013775211654,
  0.01596366921062392,
 -0.01702805985810750,
 -0.00684515010023051).
```

The interval checker proves that every one is nonzero. Thus the star
span has dimension exactly four, and (3) is unique up to scale, with
these exact algebraic cofactors as a possible choice of `lambda`.
In particular its signs are exactly `(-,+,+,-,-)` in the original
vertex order. The decimals above only identify the certified intervals.

The checker also proves that the row reduction using the first column
as pivot has exactly the degree sequence

```text
(5,4), (4,4), (3,4), (3,3), (2,3), (2,2), (1,2), (1,1).
```

Every pivot interval excludes zero. Each discarded leading coefficient
is exactly zero by the leading-vector argument in Section 2, rather
than by a numerical tolerance. Consequently the actual algebraic root
has the balanced case of (2): every `A_a,B_a` is affine.

Write

```text
h = U_22 - i U_21,       g = -U_12 + i U_11.
```

Then, exactly,

```text
Q_a = h A_a + g B_a,
Im(conjugate(h) g) = 1,
deg h = 4,       deg g = 3.                                (4)
```

The degree upper bounds follow from the seven specified row operations;
the intervals prove that the degree-four and degree-three real leading
coefficients are nonzero. Both leading coefficients are real. These
statements give an exact algebraic representation of the certified
profile; they do not give a minimal polynomial for its coefficient field.

Let `c=1/lc(h)`, which is nonzero. Monicity of every `Q_a` and
`deg g=3` force every `A_a` to have the same slope `c`. Equation (1)
then forces the slopes of the `B_a`. There are exact real algebraic
numbers `r_a,v_a,b` for which

```text
A_a = c(t-r_a),
B_a = (b+Y_a/c)t+v_a.                                     (5)
```

Thus the ten quadratic norms can be recovered from the much smaller
pencil by

```text
n_ab = (A_a B_b-B_a A_b)/(Y_b-Y_a).                        (6)
```

## 4. Consequences for the unique relation and rationality

Since all quintics are monic, (3) gives `sum lambda_a=0`.
Taking the bracket with any `Q_a` and using (1) gives

```text
sum_b lambda_b (Y_b-Y_a)n_ab = 0.                          (7)
```

Comparing leading coefficients in (7) then yields
`sum lambda_a Y_a=0`. In the balanced representation (5), the complete
constant-relation conditions are

```text
sum lambda_a = sum lambda_a Y_a
             = sum lambda_a r_a = sum lambda_a v_a = 0.    (8)
```

For example, the coefficient of `t` and the constant coefficient in
(7) also give, for every `a`,

```text
sum_(b!=a) lambda_b (Y_b-Y_a) Re(u_ab) = 0,
sum_(b!=a) lambda_b (Y_b-Y_a) |u_ab|^2 = 0.
```

The projective relation vector `[lambda_0:...:lambda_4]` is permuted by
reanchoring: a common invertible affine substitution in (3) preserves
its kernel. Hence its unordered projective class is an exact `S_5`
invariant of this profile. Rational root coordinates would make this
projective vector rational, but the converse is not established.

If all original roots were in `Q(i)`, the construction above would use
`K=Q`; the reductions and all their nonzero divisions would stay rational.
Then `h,g` would be in `Q(i)[t]`, all `A_a,B_a` in `Q[t]`, and every
pencil minor would have negative-square discriminant:

```text
disc(A_a B_b-B_a A_b)
  = -4 (Y_b-Y_a)^2 (Im u_ab)^2.                            (9)
```

This provides a more structured rational construction target: a real
rational affine pencil, a rational complex pair (4), and the ten oriented
Gaussian rational contact points at which the corresponding two
quintics vanish. The contact requirement is essential; prescribing only
the ten quadratic minors is insufficient. The five remaining private
linear factors must also be present and all root/conjugate separations
must hold.

## 5. A denominator-inclusive norm condition for a rational pencil

Apply a constant diagonal row operation and a shear to (5). This preserves
all minors and puts the affine pencil in the form

```text
A'_a=t+alpha_a,       B'_a=Y_a t+beta_a,
d_a=beta_a-Y_a alpha_a.
```

Explicitly, `alpha_a=-r_a`, `beta_a=c v_a-c b alpha_a`; the corresponding
complex frame is `h'=c h+b g`, `g'=g/c`. All these operations remain
rational if the input is rational.

For `delta=Y_b-Y_a` and `H=alpha_b-alpha_a`, direct expansion gives

```text
M_ab = A'_a B'_b-B'_a A'_b
     = delta(t+alpha_a)(t+alpha_b)
       +d_b(t+alpha_a)-d_a(t+alpha_b),
disc M_ab = (delta H+d_a+d_b)^2-4d_a d_b.                 (10)
```

Negative discriminants force all `d_a` to be nonzero with one common
sign. When that sign is positive, (10) gives the strict inequalities

```text
(sqrt(d_a)-sqrt(d_b))^2
   < delta(alpha_a-alpha_b)
   < (sqrt(d_a)+sqrt(d_b))^2.
```

In particular `alpha_a` decreases as `Y_a` increases. The direction
reverses when all `d_a` are negative.

For the certified solution the interval checker proves `d_a>0` for all
five labels and

```text
alpha_0 > alpha_3 > alpha_1 > alpha_2 > alpha_4.
```

The corresponding approximate values of `d_a`, in vertex order, are
`(0.04468125, 0.00230020, 0.07218859, 0.00339775, 0.84119327)`.
The strict signs and ordering are exact interval conclusions.

For a rational Gaussian profile, write `disc M_ab=-4s_ab^2` with
`s_ab` nonzero rational. Equation (10) then gives the exact norm identity

```text
d_a d_b = ((delta H+d_a+d_b)/2)^2+s_ab^2.                 (11)
```

Consequently all five `d_a` lie in a single class in
`Q^*/Norm(Q(i)^*)`. Equivalently, the parity of `v_p(d_a)` is independent
of `a` at every rational prime `p=3 mod 4`; negative valuations are
included. More explicitly, set `kappa=d_0`, `w_0=1`, and for `a!=0`
put

```text
w_a = (((Y_a-Y_0)(alpha_a-alpha_0)+d_0+d_a)/2+i s_0a)/d_0.
```

Then `d_a=kappa Norm(w_a)`. Evaluating (6) at `t=-alpha_a` also gives

```text
Norm(alpha_a+u_ab)
 = -d_a (alpha_b-alpha_a)/(Y_b-Y_a).                     (12)
```

These constraints retain rational denominators. They do not by
themselves imply an obstruction: for rational nonzero `kappa`, taking
`alpha_a=-kappa Y_a`, `beta_a=kappa` gives
`d_a=kappa(1+Y_a^2)` and `M_ab=(Y_b-Y_a)(t^2+kappa^2)`.
Every discriminant is a negative rational square, but all ten norm
factors coincide. Thus distinct contact factors, the rank-four
condition, and the common frame must remain in the rationality problem.

The pairwise norm identities must not be turned into an unproved global
Gaussian Gram representation. Define the real symmetric matrix

```text
R_ab = ((Y_b-Y_a)(alpha_b-alpha_a)+d_a+d_b)/2,
R_aa = d_a.
```

Since `beta_a=d_a+Y_a alpha_a`, it satisfies

```text
R = (beta 1^T+1 beta^T-Y alpha^T-alpha Y^T)/2,
```

so its rank is at most four and the relation vector lies in its kernel.
Nevertheless it need not be positive semidefinite. For the actual
certified solution the checker proves

```text
det R_{ {0,1,2},{0,1,2} } < -10^(-9)
```

(the approximate value is `-1.50970278e-9`). Thus even this real witness
has no representation `R_ab=Re(conjugate(w_a)w_b)` by one coherent
set of complex numbers. The individual norm identities (11) alone
cannot supply such a representation or its extra compatibility equations.

The coefficient field, and even whether this particular isolated root
is rational, remain undetermined. The representation does not rule out
other rational solutions or more general integer profiles.

For a fixed rational pencil and oriented contacts, the later
[eight-contact reconstruction theorem](five_star_contact_frame_reconstruction.md)
shows that an existing genuine frame and the five private roots must
also be rational. Finding compatible rational input remains necessary.
The [integer five-star lift](four_row_integer_five_star_lift.md) gives
the exact numerical counterpart of the bracket identities, including
all correction factors; it does not supply polynomial degree reduction.

Run both exact checks with

```sh
python3 docs/check_four_row_real_polynomial_certificate.py
python3 docs/check_four_row_star_pencil.py
```
