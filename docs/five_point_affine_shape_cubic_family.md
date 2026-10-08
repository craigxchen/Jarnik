# Five primitive endpoint points whose radius is cubic in affine-shape height

There is an explicit five-point family for which

```text
R ~ (4/703125) H_shape^3,
```

where `H_shape` is the maximum absolute primitive joint Plucker
coordinate: all ten triangle determinants divided by their integer
gcd. Thus no universal bound `R<=constant H_shape^2` holds, even for
five primitive centered-circle lattice points on a fixed endpoint
arc. In fact this family lies on `5 sqrt(R)` arcs, with normalized
arc constant tending to `2sqrt(5)`.

The constant does not tend to zero. This example does not refute an
inequality restricted to sufficiently small endpoint constants, does
not provide a growing number of rows, and does not contradict the
uniform endpoint-count conjecture.

It gives a lower bound of three on any universal affine-shape
reconstruction exponent, to compare with the proved exponent five
for five points in
[affine_shape_radius_divisibility.md](affine_shape_radius_divisibility.md).
The exact joint relation lattice uses the same primitive Plucker
height as in
[joint_affine_relation_lattice.md](joint_affine_relation_lattice.md).

## 1. Five actual Gaussian points

Let `T` run through the positive multiples of `600`. Use the five
subsets of `{1,...,6}`

```text
S_0={4,6},       S_1={1,3,6},       S_2={1,4,5},
S_3={2,3,5},     S_4={1,2,3,4}.
```

Every subset has sum ten. Define

```text
w_j(T)=(-1)^|S_j| product_(a in S_j)(a+iT)
                    product_(a notin S_j)(a-iT),
z_j(T)=w_j(T)/720.                                   (1)
```

Since `T` is divisible by every `a=1,...,6`, every factor `a+iT`
or `a-iT` is divisible by the ordinary integer `a`. Hence (1) is
Gaussian integral. All five points have common radius

```text
R(T)=sqrt(product_(a=1)^6(T^2+a^2))/720
    ~ T^6/720.                                      (2)
```

These points are already Gaussian-primitive. To prove this, put
`T=60n`, with `10|n`, and `b_a=60/a`, so the six coefficients are

```text
{60,30,20,15,12,10}.
```

After division by `720`, the factors are `1+i b_a n` and
`1-i b_a n`, up to the displayed row signs. Every `a` occurs both
inside and outside the five subsets. If a Gaussian prime divided
all five rows, it would therefore have to divide at least two
distinct factors from this list: one factor alone is absent from
some row.

For two factors with distinct coefficients, eliminating `n` shows
that the prime divides `b_a-b_c` or `b_a+b_c`; for opposite factors
with the same coefficient, it divides two. The prime divisors of
these finitely many nonzero integers belong to
`{2,3,5,7,11}`. The primes `3,7,11` are inert in `Z[i]` and cannot
divide a Gaussian integer of real coordinate one. Only primes over
two or five remain. But `10|n` makes every factor congruent to one
at those primes. This contradiction proves

```text
gcd_G(z_0,...,z_4)=1.                                (3)
```

In fact the same elimination proves the stronger statement that the
six block norms `Norm(1+i b_a n)` are pairwise coprime and each
block is coprime to its conjugate. Any common rational norm prime
would divide a choice of signed factors from the two blocks, which
the argument just excluded. Thus every prime valuation column of
the five-point tuple takes only two levels, zero and its full prime
exponent. By
[the monomial valuation-content floor](monomial_basis_valuation_content_floor.md),
the squared radius in (2) divides the primitive squared radius
after every nonsingular integer affine monomial row change. There
is no hidden further reduction of the radius in (2), even throughout
that entire class of row changes.

## 2. Exact triangle determinants and their full gcd

The real and imaginary coefficients of `w_j` are conveniently given
by

```text
Re w_j = -T^6-45T^4-c_j T^2+e_j,
Im w_j = -T^5-d_j T^3-f_j T,

j    c_j    e_j     d_j      f_j
0    404     720      85      1164
1    428    -720      61       396
2    464    -720      25       324
3    484    -720       5      -276
4    524     720     -35     -1236.
```

For sorted triples, write
`D_ijk(w)=det(w_j-w_i,w_k-w_i)`. Direct exact expansion gives

```text
D_ijk(w)=-960 T q_ijk(T),

ijk   q_ijk(T)
012   27T^2+108
013   56T^2+1008
014   144T^2+3600
023   50T^2+900
024   225T^2+3600
034   200T^2+3600
123   21T^2
124   108T^2+108
134   112T^2+1008
234   25T^2+900.                                     (4)
```

All ten determinants are nonzero for positive `T`, so the five
points are distinct and every three are noncollinear.

For `6|T`, the gcd of the ten integers `q_ijk(T)` is exactly `36`.
All are divisible by `36`. Write `T=6s` and divide them by `36`.
Four times the entry `012` minus entry `124` is nine, while entry
`024` is `225s^2+100`, which is not divisible by three. Thus the
remaining gcd is one.

Since division of every point by `720` divides every triangle
determinant by `720^2`, (4) proves the exact difference-lattice index

```text
g(T)=gcd_(i<j<k) |D_ijk(z)|
    =34560T/720^2=T/15.                              (5)
```

The largest entry in (4) for `T>=1` is `q_024`. Therefore

```text
H_shape(T)=(225T^2+3600)/36
          =25T^2/4+100.                              (6)
```

The primitive joint height, not separate projective heights of
individual affine relations, is what appears in (6).

### The normalized affine conic is explicit

Set `x=T^2` and send `z_0,z_1,z_2` affinely to
`(0,0),(1,0),(0,1)`. The remaining two points become

```text
(u,v)=(-50(x+18)/(27(x+4)), 56(x+18)/(27(x+4))),
(s,t)=(-25(x+16)/(3(x+4)), 16(x+25)/(3(x+4))).
```

The unique positive conic has equation

```text
A X^2+2B XY+C Y^2-A X-C Y=0,
A=4(x+25)(x+36),
B=10(x^2+43x+360),
C=25(x+9)(x+16),
AC-B^2=(180T)^2,
A+C-2B=9x(x+1).
```

The actual Gram matrix of `z_1-z_0,z_2-z_0` is
`((x+4)/3600) ((A,B),(B,C))`. Its circumradius formula recovers
exactly `R^2=product_(a=1)^6(x+a^2)/720^2`. Thus the counterexample
also has a direct conic reconstruction with its metric scale and
Gaussian integrality both retained.

## 3. Radius growth and the endpoint arc

Combining (2) and (6) gives

```text
R(T)/H_shape(T)^2 ~ T^2/28125 -> infinity,
R(T)/H_shape(T)^3 -> 4/703125.                        (7)
```

Thus every proposed universal reconstruction exponent smaller than
three is contradicted by this family.

For `T>=600`, both coordinates of every `z_j` are negative. Their
real coordinates strictly decrease with `j`, and their imaginary
coordinates strictly increase. Their arguments are consequently
ordered on one arc shorter than a quadrant, with endpoints `z_0,z_4`.
The endpoint difference is exactly

```text
z_4-z_0=(-T^2/6)+i(T^3/6+10T/3),
d(T)=|z_4-z_0|=(T/6)sqrt(T^4+41T^2+400).             (8)
```

The short-arc length and its normalized value are therefore

```text
ell(T)=2R(T) arcsin(d(T)/(2R(T))),
C(T)=ell(T)/sqrt(R(T)) -> 2sqrt(5).                  (9)
```

For an explicit uniform bound, `T>=600` gives

```text
d(T)/sqrt(R(T)) <= 2sqrt(5) sqrt(1001/1000),
d(T)/(2R(T)) <= 1/1000.
```

Use `arcsin u<=u/sqrt(1-u^2)` to conclude `C(T)<5`.
This verifies the claimed endpoint condition without a numerical
angle test. The positive limiting constant is retained: no
arbitrarily small-constant family is asserted.

## 4. The forward polynomial and why reversal matters

Before reversal, the polynomial points are

```text
v_j(t)=(-1)^|S_j| product_(a in S_j)(at+i)
                    product_(a notin S_j)(at-i),
w_j(T)=T^6 v_j(1/T).
```

The common subset sum makes all five `v_j` agree through order two
at `t=0`. Their triangle determinants are exactly
`-960 t^9 q_ijk(1/t)t^2`; the common polynomial divisor has order
nine, not merely the order-six contribution from two differences.
Reversal turns these into (4), and turns the small-parameter
coalescence into the endpoint family at large integer `T`.

For example the forward `123` determinant is exactly
`-20160t^9`, certifying that no larger common power is present.
The complete expansions, not just this heuristic about orders,
are checked in the accompanying exact verifier.

## Verification

Run `python3 docs/check_five_point_affine_shape_cubic_family.py`.
It checks the two polynomial presentations, every coefficient in
(4), the exact ninth-order forward divisor, equal norms, primitive
Gaussian gcd, the full triangle gcd and joint height, and the finite
endpoint-arc inequalities by integer arithmetic. The unbounded
conclusions in (7) follow from the explicit formulas.

## 5. Exact signed-shape obstruction to a shrinking constant

Keep the five labelled subsets from Section 1, but replace the six
coefficients `1,...,6` by arbitrary real coefficients having the same
five subset sums.  Solving those four linear equations gives exactly

```text
(a_1,...,a_6)=(u,2u,v,u+v,2u+v,3u+v).                (10)
```

Thus this parameterization includes every signed rational shape, not
only the positive chamber.  Put

```text
P=product_r a_r=2u^2 v(u+v)(2u+v)(3u+v),
s_r=sqrt(T^2+a_r^2),
w_j=(-1)^|S_j| product_(r in S_j)(a_r+iT)
                    product_(r notin S_j)(a_r-iT).
```

Assume for now that `P!=0`.  Direct cancellation of the common oriented
factors gives the two exact identities

```text
w_1-w_4=-4u(u+v)(3u+v)
             (a_1+iT)(a_3+iT)(a_5-iT),
w_2-w_4=-4uv(2u+v)
             (a_1+iT)(a_4+iT)(a_6-iT).               (11)
```

The elementary reason is that the three changed coefficients obey
`a_2+a_4=a_6` and `a_2+a_3=a_5`; the real part of the corresponding
three-factor expression is constant.  Taking absolute values in (11)
and using `|w_j|=s_1...s_6` gives

```text
|w_1-w_4| |w_2-w_4| / |w_0|
    =8|P| s_1/s_2
    =8|P| sqrt((T^2+u^2)/(T^2+4u^2)).                 (12)
```

Equivalently, the squared polynomial identity is

```text
N(w_1-w_4) N(w_2-w_4)(T^2+4u^2)
 =64P^2 N(w_0)(T^2+u^2).                              (13)
```

If `z_j=w_j/P` are Gaussian integral, their radius is
`R=|w_0|/|P|`, and (12) becomes

```text
|z_1-z_4| |z_2-z_4| / R
 =8 sqrt((T^2+u^2)/(T^2+4u^2)) >=4.                  (14)
```

Every arc containing the five points has length at least both chords.
Therefore

```text
arc length / sqrt(R) >=2.                             (15)
```

This is an all-parameter statement.  It covers arbitrary signs,
varying rational ratios, finite-ratio degenerations, and different
relative rates of `u`, `v`, and `T`; no tangent expansion or positivity
assumption is present.

For comparison, the five cubic tangent coefficients are

```text
d_0=u(13u^2+12uv+4v^2),
d_1=u(13u^2+10uv+2v^2),
d_2=u(u^2+2uv+2v^2),
d_3=5u^3,
d_4=u(u^2-6uv-2v^2).
```

They obey

```text
(d_1-d_4)(d_2-d_4)=8P,
```

so the limiting cleared constant is at least `2sqrt(2)` whenever all
six factors are tangential.  At `v/u=-3/2`, where equality occurs in
that tangent estimate, rows `1` and `2` coincide exactly.  The exact
chord estimate (15) is stronger for the present purpose because it also
handles every nontangential signed regime.

If `P=0`, cancellation of the common zero factor does not leave five
distinct rows.  The row counts after common-factor removal satisfy

| Degeneration | At most this many distinct rows |
|---|---:|
| `u=0` | 1 |
| `v=0` | 3 |
| `v=-u` | 2 |
| `v=-2u` | 2 |
| `v=-3u` | 3 |

Thus none of the exact boundary shapes evades (15) with five points.

## 6. What primitive normalization can change

If the integral tuple `z_j=w_j/P` has common Gaussian divisor `G`,
primitive reduction changes (14) to

```text
|z'_1-z'_4| |z'_2-z'_4| / R'
    >=4/|G|,
arc length / sqrt(R') >=2/sqrt(|G|).                  (16)
```

Consequently a shrinking primitive constant in this construction can
only come from an unbounded common Gaussian divisor after the displayed
clearing.

Here is one exact regime where that escape is excluded uniformly.
Assume `u,v,T` are integers, all six coefficients `a_r` in (10) are
nonzero, each `a_r` divides `T`, and use the literal integral clearing

```text
z_j=(-1)^|S_j| product_(r in S_j)(1+iq_r)
                    product_(r notin S_j)(1-iq_r),
q_r=T/a_r.                                             (17)
```

Then, with a deliberately nonoptimal constant,

```text
N(gcd_G(z_0,...,z_4)) <= 2^6 5^2=1600.                (18)
```

For the ramified prime, `v_(1+i)(1+iq)=v_(1+i)(1-iq)=1`
when `q` is odd and is zero otherwise, so the common exponent is at
most six.  An inert odd prime cannot divide a factor of norm `1+q^2`.
For an odd split prime `p`, put `E=max_r v_p(a_r)`.  If
`v_p(T)>E`, every `q_r` is zero modulo `p` and no factor is divisible.
Otherwise only coefficient forms at the maximal valuation can be
active.  After dividing the common power from `u,v`, unequal residual
valuations leave an active set contained in `{1,2}` or the singleton
`{3}`; their signed supports do not cover all five rows.  If both
residual values are units but one of `a_4,a_5,a_6` vanishes modulo
`p`, that single form is the only maximal one for `p>=5`.  Thus the
only remaining case has all six forms as units.

In that case, mark the five-row support of every factor divisible
by one chosen Gaussian prime, using the complementary support for
`1-iq_r`.  There are 41 minimal signed covers.  For the full-rank
covers, the gcds of the coefficient minors have prime support contained
in `{2,3,5}`.  The ten rank-deficient covers force one of the six
coefficient forms to vanish; primes at which their rank drops are again
in `{2,3}`.  This treats the two-factor covers as well as the covers
having three or more factors, rather than silently omitting zero
three-by-three minors.  Since `3` is inert, only `5` remains.  At that
prime the sole nondegenerate cover is

```text
factor 2 positive, factor 4 positive, factor 5 negative.
```

Indeed the mod-five solution is `v=u`, and the six coefficient forms
are `(u,2u,u,2u,-2u,-u)`.  These are the only three tagged factors at
the chosen Gaussian prime, even if one began with a nonminimal cover.

Rows `1`, `2`, and `3` are covered only by factors `5`, `4`, and `2`,
respectively.  If their three valuations were all at least `h`, inversion
of the unit coefficients would give

```text
v=u mod 5^h,       v=-4u mod 5^h.
```

After removing the common coefficient valuation, `u` is a unit in the
only case where this cover occurs, so `h<=1`.  Applying the safe bound
to both primes above five proves (18).
Combining it with (16) leaves an absolute positive lower bound for the
primitive normalized arc constant throughout (17); explicitly it is at
least `1/sqrt(10)`.

The signed identities and the finite support enumeration are checked by
[check_five_point_signed_tangent.py](check_five_point_signed_tangent.py).
The literal argument above is subsumed by the prime-by-prime content
bound in
[six_factor_primitive_content_bound.md](six_factor_primitive_content_bound.md).
That bound removes the assumptions `a_r|T` and proves the same numerical
lower bound for every rational projective parameter triple in this
six-factor family.
