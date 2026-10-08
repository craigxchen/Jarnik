# Contact divisors with conjugate isotropic directions

At a contact where `x^2+y^2` is anisotropic, a zero of a Gaussian norm
forces both ground-field components to vanish.  This fails at a closed point
whose residue field contains `i`: one isotropic factor may vanish while
the other remains a unit.  The correct replacement is a weighted
global divisor inequality that keeps both conjugate isotropic factors.

The inequality is valid on a smooth projective curve of any genus.  For
four fixed directions, effective contact divisors of degrees `7,7,7,4`
over a field in which `-1` is nonsquare still force common pole degree at
least five.  Without the odd-degree anisotropic contribution, the
bound drops to four and is sharp at the level of divisor accounting.

## 1. The two isotropic rows of a determinant-one frame

Let `k` be a characteristic-zero field in which `-1` is not a square,
let `C/k` be a smooth projective geometrically integral curve, and put
`K=k(C)`.  Formal reality is sufficient but not necessary; the proof
uses only the quadratic constant extension `k(i)/k` and the local
anisotropy condition stated below.
Take

```text
U=((a,c),(b,d)) in SL_2(K),       ad-bc=1.
```

Over `k(i)` write

```text
z_+=a+ib,       w_+=c+id,
z_-=a-ib,       w_-=c-id.                             (1)
```

These are the two conjugate isotropic rows.  Their determinant is the
constant

```text
z_+w_- - z_-w_+ = -2i.                               (2)
```

Let `Delta` be the least common pole divisor of `a,b,c,d` and put
`delta=deg Delta`.  If `sigma` is the canonical section of
`O_C(Delta)`, then

```text
Z_+=sigma z_+,       W_+=sigma w_+,
Z_-=sigma z_-,       W_-=sigma w_-                   (3)
```

are sections after base change to `k(i)`, and (2) becomes

```text
Z_+W_- - Z_-W_+ = -2i sigma^2.                        (4)
```

Now fix a ground-field direction
`lambda=[alpha:beta] in P^1(k)` and define

```text
L_(lambda,+)=alpha Z_+ + beta W_+,
L_(lambda,-)=alpha Z_- + beta W_-,

F_lambda
 =(alpha a+beta c)^2+(alpha b+beta d)^2
 =(L_(lambda,+)L_(lambda,-))/sigma^2.                 (5)
```

Neither section in (5) is identically zero.  If one were zero, its
constant conjugate would also be zero, and the nonzero ground-field
direction `(alpha,beta)` would lie in the kernel of the invertible
matrix `U`.
Consequently each section has a zero divisor of degree `delta`, and
there is an exact divisor identity

```text
Z(L_(lambda,+))+Z(L_(lambda,-))
   =2Delta+div(F_lambda).                              (6)
```

Strictly, (6) is an identity on `C_(k(i))`, with the right side pulled
back from `C`.  Equivalently, the conjugate product is a norm section
whose zero divisor descends to the displayed divisor on `C`.  Scalar
extension preserves divisor degree, so each conjugate section has
degree `delta` and their combined degree is `2delta`; there is no extra
factor of two.

At a closed point `P`, the coefficient on the right is

```text
2 ord_P(Delta)+ord_P(F_lambda).                        (7)
```

It is nonnegative because the left side is effective.  Formula (6),
not separate treatment of `a^2` and `b^2`, is the exact pole-and-zero
accounting at arbitrary closed points.

## 2. Anisotropic contacts cost four units

Call a closed point `P` **anisotropic** if `-1` is not a square in its
residue field `k(P)`, and **split** otherwise.  Suppose

```text
r=alpha a+beta c,       s=alpha b+beta d,
k_P=ord_P(r^2+s^2)>0.                                 (8)
```

At an anisotropic point there is no leading cancellation between the
two squares.  With the convention that the order of the zero function
is infinity,

```text
k_P=2 min(ord_P(r),ord_P(s))=2e_P,
e_P>=1.                                               (9)
```

Choose a complementary direction `[gamma:eta]` with
`alpha eta-beta gamma!=0`, and let `(r',s')` be its two ground-field
components.  The determinant-one identity gives

```text
rs'-sr'=(alpha eta-beta gamma)(ad-bc)
       =alpha eta-beta gamma.                         (10)
```

Since both `r` and `s` vanish to order at least `e_P`, at least one of
`r',s'` has a pole of order at least `e_P`.  Therefore

```text
ord_P(Delta)>=e_P=k_P/2.                              (11)
```

Combining (7), (9), and (11), an anisotropic norm zero of order `k_P`
uses at least

```text
2ord_P(Delta)+k_P >= 2k_P.                            (12)
```

In particular a reduced prescribed contact, which only asks for
`k_P>=1`, actually has `k_P>=2` and contributes at least four to the
left side of (6).  At a split point a reduced contact can be paid by a
simple zero of one isotropic factor and contributes only one.  This is
the exact source of the weights four and one below.

## 3. The weighted aggregate inequality

Let `Lambda` be a finite set of `r_0` ground-field directions.  For each
`lambda in Lambda`, let `D_lambda` be an effective prescribed contact
divisor such that

```text
D_lambda <= Z(F_lambda).                              (13)
```

Assume first that the `D_lambda` are reduced.  Let `A_lambda` be the
part of `D_lambda` supported at
anisotropic closed points, and put

```text
D=sum_lambda deg D_lambda,
A=sum_lambda deg A_lambda.                            (14)
```

Sum (6) over all directions, but retain on its left only zeros lying
over the prescribed contacts.  Every split contact contributes at
least one by (7).  Every anisotropic contact contributes at least four
by (12).  Since the two section-zero divisors for each direction have
total degree `2delta`, this proves

```text
D+3A <= 2r_0 delta.                                   (15)
```

Anisotropic contacts for distinct directions are automatically
disjoint.  Indeed, if the norms for two distinct ground-field
directions vanished at the same anisotropic point, (9) would make both
ground-field component pairs vanish there.  The two directions form an
invertible constant change of basis, so all four entries of `U` would
vanish there, contradicting `det(U)=1`.  Equation (11) now puts each
anisotropic contact at a distinct point of `Delta`, and hence

```text
A<=delta.                                             (15a)
```

For rational contacts this second inequality recovers the sharper
pointwise conclusion: if all twenty-five contacts are rational, then
`A=25` and `delta>=25`.  The aggregate inequality (15) is useful when
most contacts are split and do not themselves force poles.

There is also a multiplicity version.  If a contact `P` occurs in
`D_lambda` with multiplicity `m`, its split weight is `m`.  At an
anisotropic point, (9) forces
`k_P>=2ceil(m/2)`, and (11)--(12) give weight at least

```text
4ceil(m/2).                                           (16)
```

Thus, if `m_(lambda,P)` denotes the coefficient of `P` in
`D_lambda`, the multiplicity-weighted sum

```text
W = sum_(lambda,P split) m_(lambda,P)deg(P)
    +sum_(lambda,P anisotropic)
       4ceil(m_(lambda,P)/2)deg(P)
  <=2r_0delta.                                         (16a)
```

This remains valid with overlapping supports because contacts are
counted with their direction incidence.  Since anisotropic supports
are automatically disjoint across directions, the multiplicity form
of the separate pole bound is also unconditional:

```text
sum_(lambda,P anisotropic)
  ceil(m_(lambda,P)/2)deg(P) <= delta.                 (16b)
```

The inequality has no genus term.  Genus can impose separate gonality
or Riemann--Roch restrictions on the two maps
`[z_+:w_+]` and `[z_-:w_-]`, but it does not change the local weights in
(15).  This makes (15) portable to rational, elliptic, and higher-genus
fixed sources.

## 4. The `7+7+7+4` profile when `-1` is nonsquare

Take the four directions

```text
z,       w,       z+w,       -3z+2w.                  (17)
```

Suppose their effective contact divisors, with no reducedness or
disjointness assumption, have degrees

```text
7,       7,       7,       4.                         (18)
```

A split closed point has residue field containing `k(i)`.  Since
`k(i)/k` is quadratic, its degree over `k` is therefore even.  In each
of the first three contact divisors, the total degree is odd while the
split part has even degree.  Some anisotropic point `P` must therefore
occur with both odd multiplicity `m` and odd residue degree.  At this
incidence, the multiplicity weight in (16a) exceeds its ordinary
divisor-degree contribution by

```text
[4ceil(m/2)-m]deg(P)=(m+2)deg(P)>=3.                  (19)
```

The three odd-degree directions supply three such excesses.  Applying
(16a) with `r_0=4` gives

```text
34<=W<=8delta,
delta>=5.                                             (20)
```

Thus the degree-four short-frame threshold is still excluded even when
the contacts are arbitrary closed points, provided the four contact
divisors are defined over a field with `-1` nonsquare and have the odd
degree pattern (18).  On an elliptic source, the projective frame map
`[1:a:b:c:d]` has pullback degree `delta`, so (20) yields the same
height gap as the rational-contact argument.

This is a pointwise function-field degree obstruction: the frame and
the contact divisors may vary from one instance to another.  Turning
it into a lower bound for heights after evaluating a sequence of frame
maps requires either a fixed map or independent control of the moving
coefficients, as in the
[bounded-degree moving-frame height lemma](moving_bounded_degree_frame_contact_height.md).

This is a degree-one profile obstruction.  More explicitly, on a fixed
elliptic curve over `Q`, effective contact divisors over `Q` of degrees
`7,7,7,4` rule out a rational determinant-one frame map of degree at
most four, even if the divisors are nonreduced, overlap, or contain no
rational contacts.  If all four degrees are instead scaled by an
integer `q`, odd `q` gives only three parity excesses and hence
`W>=25q+9`; this beats `8delta<=32q` only for `q=1`.  Even `q` gives no
parity excess.  In the reduced case, writing `A` for the total
anisotropic degree gives the more informative bound

```text
25q+3A <= 8delta.                                      (20a)
```

To contradict the correspondingly scaled threshold `delta<=4q`, one
needs

```text
A>7q/3.                                                (20b)
```

For odd `q`, parity alone says only `A>=3`, in agreement with the
three-excess calculation above.  Thus the argument does not by itself
exclude scaled degree-four profiles without additional anisotropic
contact information.

This limitation is realized over `Q` already at `q=2`: the
[explicit rational frame](scaled_split_contact_profile_q2.md) has
common pole degree seven and reduced, disjoint, all-split contact
divisors of degrees `14,14,14,8`. It satisfies the four norm-contact
conditions, without asserting the full source or lattice-arc profile.

This conclusion uses divisors over `k` of the degrees in (18), not
merely twenty-five geometric points after base change.  If the contact
data are specified only over an algebraically closed field, the parity
contribution in (19) disappears.

## 5. Elementary examples and sharpness

The basic failure of the pointwise anisotropic argument is

```text
C=P^1,       U(t)=((t,-1),(1,0)).                     (21)
```

Here `Delta=[infinity]`, while the norm of the first column is
`t^2+1`.  Its degree-two closed zero contains the two conjugate
isotropic geometric directions: at `t=i` one of `t+i,t-i` vanishes,
and at `t=-i` the other vanishes.  No frame entry has a pole there.
In (6) these two simple isotropic zeros use the full degree two of the
two conjugate sections, with no contact pole contribution.

A four-direction version shows that the coefficient `2r_0` in (15)
cannot be reduced by determinant one alone.  Let `f in k(P^1)` be a
polynomial of degree `m` and take

```text
U_f=((f,f-1),(1,1)),       det(U_f)=1.                (22)
```

Its common pole divisor is `m[infinity]`.  The four norm functions in
(17) are

```text
f^2+1,
(f-1)^2+1,
(2f-1)^2+4,
(f+2)^2+1.                                            (23)
```

Their zero divisors are pairwise disjoint for the displayed distinct
target values, each has degree `2m`, and every zero is isotropic; for a
general `f` these complete zero divisors need not be reduced.  Thus
the total contact degree over all four complete zero divisors is

```text
8m=2r_0 deg Delta,                                    (24)
```

so (15) is an equality with `A=0`.  For `f=t^m`, all eight target values
are nonzero, so the fibers are reduced.  In particular, for `m=4` the
common pole degree is four and the complete split contact degree is
thirty-two.  Over an algebraic closure one may select disjoint reduced
subsets of sizes
`7,7,7,4` from those fibers without forcing another pole.  Those
subsets do not descend to the odd-degree divisors in (18): conjugate
isotropic points descend in even-degree packets.  This example both
verifies the accounting and shows why the anisotropic parity term is
essential.

The weighted identity (15) is therefore sharp for arbitrary split
contacts.  It extends the rational-contact pole argument exactly as
far as determinant one and constant conjugation allow; any stronger
claim for even-degree or algebraically closed contact profiles needs
additional geometry beyond the frame identity.
