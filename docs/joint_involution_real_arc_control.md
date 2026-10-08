# Choosing a joint involution that preserves the original arc

For six distinct real circle points in a proper arc, two opposite ordered
labelings suffice to choose an arc-preserving version of the joint Segre
involution, provided its relevant matching difference is nonzero. That
nonvanishing is automatic for the increasingly accurate full-profile
arithmetic inputs. The new six points lie inside the subarc spanned by
the first four points in the chosen ordering. There is no angular loss
and no lower bound on the spacing of the old points is needed.

The radius is a separate arithmetic quantity. The subsequent sixteen-cut
argument proves that this particular operation increases the least integral
radius even without retaining the original anchor. Its normalized output
arc length becomes unbounded for every labeling of the fixed-three pullback.

## 1. The exact pole in the inverse normalization

Choose an increasing real projective coordinate `x` on the original
proper arc, and translate it so the first point is `x_1=0`. Write
`d=x_2`, `e=x_3`, and `L=x_6`, so `0<d<e<L`. Normalizing the
first three points to `infinity,0,1` gives

```text
t=chi(x)=e(x-d)/((e-d)x),
tau=e/(e-d),        kappa=ed/(e-d),
x=chi^(-1)(t)=kappa/(tau-t).                             (1)
```

The point `t=tau` is the pole of this affine inverse. For the remaining
normalized old coordinates one has

```text
1<a=chi(x_4)<b=chi(x_5)<c=chi(x_6)<tau,
kappa=L(tau-c).                                         (2)
```

Consequently a general proposed normalized output point `t` satisfies

```text
|x(t)|=L/rho(t),
rho(t)=|tau-t|/(tau-c).                                  (3)
```

This is the precise separation parameter needed for a fixed labeling.
If all new finite normalized points satisfy `rho>=rho_0>0`, their
affine coordinate sizes are at most `L/rho_0`. For example, for an
original angular span `Theta<pi`, take
`x=tan((theta-theta_1)/2)` and `L=tan(Theta/2)`. A two-sided
bound then places the new angles in an interval of length at most

```text
4 arctan(tan(Theta/2)/rho_0).                             (4)
```

If all the new coordinates have `t<tau`, their inverse coordinates
are positive, and the factor four in (4) improves to two.
No rational-height estimate alone gives this real pole separation.

## 2. An ordered choice avoids the pole exactly

Use the formulas from
[joint_segre_involution.md](joint_segre_involution.md):

```text
F=c(a+1-b)-a,      G=c(a+b-1)-a,
H=c(b+1-a)+a-2b,   U=c(b+1-a)-a,
V=c(a+b-1)-a(2b-1),
(a',b',c')=(aH/V,bF/V,cH/F).                             (5)
```

In the ordered chart `1<a<b<c`, all five original matching
coordinates are positive. Therefore `H,U,V,G` are positive by their
expressions as positive linear combinations of those matchings.
Assume `F<0`. The exact collision identities imply

```text
c'<b'<0,            1<a'<a.                              (6)
```

Indeed `b',c'<0`, and
`c'-b'=(c-b)UG/(FV)<0`. Also
`a'-1=-(a-1)F/V>0`, while
`a'-a=a(H-V)/V=-2aA/V<0`.

The inverse normalization is increasing on `t<tau`. Negative `t`
corresponds precisely to points strictly between old rows `1` and
`2`, and `1<t<a` corresponds to the interval between old rows `3`
and `4`. Thus the exact output order on the original arc is

```text
old1 < new6 < new5 < old2 < old3 < new4 < old4.            (7)
```

All six output points lie in the original first-to-fourth subarc.
In particular they lie in the original arc, regardless of how uneven
the old spacings are. This conclusion is projective-order based, so
it applies to any proper circle arc after choosing an affine parameter
whose pole lies outside it. The small-angle chart used for (4) is not
needed for (7).

If the initially ordered chart has `F>0`, reverse all six labels and
normalize the first three in the reversed order. The positive ratio of
the matching coordinates `B/C` becomes `C/B`: reversal takes the raw
matching products `B,C` to `-C,-B`, and normalization changes both by
one common scalar. Therefore the reversed ordered chart has `F<0`.
One of the two opposite orderings always satisfies (6), unless `F=0`.

In the full arithmetic profile, the ratio `B/C` has rational height
`12w+o(w)`, so it cannot equal one. Hence the exceptional case `F=0`
does not occur for sufficiently accurate extracted inputs. The two
choices are conjugates of the joint map by label reversal; no claim
that there are only two conjugates in its full symmetry orbit is needed.

This proves an actual arc-preserving choice, not a conjectural choice
among a large family of transformations. It does not assert a universal
strict contraction factor for the angular diameter: the first three
fixed points alone can already span nearly the entire original arc.

## 3. A fixed labeling really can hit the pole

Pole separation would be necessary without the ordering choice. Take
the fixed normalized six-tuple

```text
(infinity,0,1,2,5/2,5).
```

Here `F=1/2`, and the new third moving coordinate is `c'=45`.
For any positive rational `epsilon`, realize these six points in
the real coordinate

```text
x=40 epsilon/(45-t).
```

Their actual coordinates are

```text
0, (8/9)epsilon, (10/11)epsilon,
(40/43)epsilon, (16/17)epsilon, epsilon.
```

As circle coordinates `x=tan(theta/2)`, they lie in an arc of angular
length `2 arctan(epsilon)`. But the same fixed labeling sends one
output to `x=infinity`, hence to the antipodal circle point. Thus
there is no bound by a fixed multiple of the original angular span
for this labeling as `epsilon->0`.

Reversing the labels gives the ordered normalized chart
`(a,b,c)=(9/4,3,6)`, with `F=-3/4`. Its three lifted new coordinates
are, in the same original `x` parameter,

```text
(250/269)epsilon, (312/331)epsilon, (160/161)epsilon,
```

and the retained three are
`epsilon,(16/17)epsilon,(40/43)epsilon`. All are inside the original
arc, as (7) requires. This is a real-geometric example; it does not
assert that these specializations have the full arithmetic profile
or a bounded endpoint constant at their least integer radius.

## 4. Actual radius inflation

Let the original angular diameter be `Theta`, with radius `R` and
normalized arc constant `C=Theta sqrt(R)`. The ordered choice above
gives `Theta_out<=Theta`. At the least integral output radius,

```text
C_out=Theta_out sqrt(R_out)
       <=C sqrt(R_out/R).                                (8)
```

Thus an additional upper bound `R_out<=R` would preserve the endpoint
constant. Arc preservation alone does not provide that upper bound.
It does strengthen the known radius lower bound. The common-phase and
integer Bezout argument for chord matchings gives, for any integral
realization of these actual output angles,

```text
h_out<=3log(R_out Theta_out).
```

For the finite radius comparisons choose `w=W/64`, where `W=2log R`
is the original inherited norm scale; the centrally truncated block
weights satisfy the stated `(1+-eta)w` bounds with this choice.
Thus `Theta_out<=Theta_in<=C exp(-16w)` exactly.
Combining this with the actual lower bound `h_out>=42w-o(w)` gives

```text
log R_out>=30w-o(w).                                     (9)
```

More precisely, using `beta=2s+log T` and the error `eta` in the joint
involution note, the lower bound is
`30w-114eta w-80beta-(16/3)log2-log C`. There is no assumption that
the output already satisfies an endpoint bound relative to `R_out`.
The weaker `28w` floor applies to arbitrary endpoint-scale circle
realizations of the same projective moduli, whose actual angular arcs
need not be the one preserved here.

The inherited input scale is `log R=32w+o(w)`. The preliminary floor (9)
is now strengthened by the actual common Gaussian matching content.
The [sixteen-cut argument](joint_involution_standalone_inflation.md) gives

```text
log R_out^(6)>=(98/3)w-o(w),
R_out^(6)>=R^(49/48-o(1)),
Theta_out sqrt(R_out^(6))>=exp(w/3-o(w)).                 (9a)
```

Thus the standalone six-point operation necessarily inflates the radius
and loses a bounded normalized endpoint constant. Retaining the original
anchor cannot reduce the necessary radius, so the same lower bounds apply
to `R_out^(7)`.

The last bound in (9a) in fact holds for every labeling of the fixed-three
pullback. Its proof combines `log(R_out Theta_out)>=(50/3)w-o(w)` with
the unchanged three directions' `Theta_out>=exp(-16w+o(w))`; it does not
require the arc-preserving upper bound. These conclusions do not address
a subsequent arbitrary fractional-linear change or another joint map.
The uniform-bound goal remains unproved.

## 5. Two retained-anchor denominator consequences

There is a further exact consequence in the original centered Gaussian
circle model. It is important here to distinguish the six output points
alone from those six points together with the original anchor `0`.
The joint map retains three of the six nonanchor directions.

Let `R_3` be the least integral radius for just this retained triple.
At a split prime this radius uses half the range of their three original
allocation exponents. Exactly forty-eight of the sixty-four core cuts
distinguish these three rows. The central correction can increase the
range by at most `2r_p`. Consequently

```text
24(1-eta)w<=log R_3<=24(1+eta)w+s.                         (10)
```

Divide the three old Gaussian points by their common Gaussian gcd to
obtain a primitive triple on this least circle. Any centered integral
circle realization retaining their relative angular directions is this
triple multiplied by one Gaussian integer `lambda`: the multiplier is
in `Q(i)`, and an integral Bezout combination of the primitive triple
shows that it lies in `Z[i]`.

In the complex coordinates of that primitive triple circle, write the
three proposed new points as reduced Gaussian fractions `a_j/b_j`.
They are rational Gaussian numbers because their half-angle ratios are
rational. The least extension multiplier for the six output directions
alone is exactly

```text
lambda_6=lcm_G(b_4,b_5,b_6),       R_out^(6)=R_3 |lambda_6|.
```

Equations (9a)--(10) therefore imply

```text
log|lambda_6|>=(26/3)w-o(w),
max_j log|b_j|>=(26/9)w-o(w).                            (11)
```

For the second assertion use `|lambda_6|<=|b_4 b_5 b_6|`.

If the original anchor `0` is also retained, the fixed set has four
directions. A core cut distinguishes this set precisely when it meets
one of the three retained nonanchor rows. There are fifty-six such cuts,
so its least radius `R_{0,3}` instead satisfies

```text
28(1-eta)w<=log R_{0,3}<=28(1+eta)w+s.                    (12)
```

Normalize this four-point set to a primitive Gaussian tuple, and write
the same three moving directions on its circle with reduced denominators
`btilde_4,btilde_5,btilde_6`. The identical Bezout argument gives

```text
lambda_7=lcm_G(btilde_4,btilde_5,btilde_6),
R_out^(7)=R_{0,3}|lambda_7|,
log|lambda_7|>=(14/3)w-o(w),
max_j log|btilde_j|>=(14/9)w-o(w).                        (13)
```

Here (13) uses (9a) for the seven-point extension, subtracting `28w+o(w)`
in (12). These are different denominator normalizations: the `(26/3)w` cost
in (11) is relative to a primitive triple, whereas the `(14/3)w` cost in (13)
is relative to the primitive fixed four-point set. Both force a substantial
extension beyond the indicated retained
circle. Neither asserts that these denominator primes are new to the full
old circle; they may come from its old conductor. These are specifically
centered Gaussian-circle statements, unlike the center-independent
matching-height inequality used for (9).

## Audit

The uniformity-audit agent independently checked the signs in (6),
the order interpretation (7), and the reversal argument. The exact
fractions in the fixed-pole example and its reversed output were checked
using the rational joint-map formulas and inverse normalization. The
conclusion concerns real angular containment, with the radius issue
retained explicitly in (8).
