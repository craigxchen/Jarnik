# Moving six-factor coefficients are excluded at normalized arc constant `1/2`

This note excludes the exact five-row signed six-factor family from the
`C<=1/2` test even when its coefficients vary. The earlier positive lower
bound already excluded normalized constants tending to zero within this
family; the refinement here reaches the specified test constant. It does
not address arbitrary circle tuples or improve the general growth bound.

Use the classified coefficient forms

```text
(a_1,...,a_6)=(u,2u,v,u+v,2u+v,3u+v),
```

and the subsets

```text
S_0={4,6}, S_1={1,3,6}, S_2={1,4,5},
S_3={2,3,5}, S_4={1,2,3,4}.
```

For an integer triple `(u,v,T)` with `gcd(u,v,T)=1`, define

```text
P=product_j a_j,
w_i=(-1)^|S_i| product_(j in S_i)(a_j+iT)
                 product_(j notin S_i)(a_j-iT).
```

Assume first `P!=0`, let `G` be the full Gaussian gcd of the five rows,
and write `|G|=sqrt(N(G))`. No distinct-absolute-value hypothesis is needed
for this representative family's content proof. The separate code
classification transfers the result to other equal-sum codes under its
stated nonzero, distinct-absolute-value hypotheses.
The later [repeated-magnitude classification](sixfactor_repeated_magnitude_classification.md)
also identifies both repeated-coefficient cases permitted by short-arc
allocation rigidity with this representative. Thus the exact equal-first-
moment six-factor exclusion now includes repeated magnitudes as well.

The existing odd-prime valuation argument in
[six_factor_primitive_content_bound.md](six_factor_primitive_content_bound.md)
gives the safe bound

```text
odd-part(|G|/|P|) <= 5.                              (1)
```

The factor `5` is the sole possible extra split-prime contribution; all
inert primes and all split primes other than `5` are charged to `P`. We only
need this deliberately non-optimized odd-prime estimate.

## The exact ramified excess

Assume `T != 0`, and put

```text
s_j=v_2(a_j),   t=v_2(T).
```

The valuation of either orientation of a factor at `1+i` is

```text
v_(1+i)(a_j plus-or-minus iT)=2 min(s_j,t)+1_(s_j=t).
```

The orientation is irrelevant at this prime, so the common gcd has

```text
Delta_2 := v_(1+i)(G)-v_(1+i)(P)
         = E-2D,
E=# {j:s_j=t},
D=sum_(s_j>t)(s_j-t).                              (2)
```

We claim `Delta_2<=1` for every primitive `(u,v,T)`. The proof is a short
parity case split.

If `t=0` and `u,v` are odd, the odd forms are `a_1,a_3,a_5` and the even
forms are `a_2,a_4,a_6`. Here `s_2=1`, while
`a_6-a_4=2u` has valuation one, so one of `a_4,a_6` has valuation one and
the other has valuation at least two. Hence `E=3` and `D>=4`. If `u` is
odd and `v` even, the same argument uses the even forms `a_2,a_3,a_5`
and `a_5-a_3=2u`. If `u` is even and `v` is odd, the four odd forms are
`a_3,a_4,a_5,a_6`, while `s_1>=1` and `s_2=s_1+1>=2`; hence `E=4` and
`D>=3`. If both `u,v` are even, then `E=0`. Thus (2) is at most one in
all `t=0` cases.

Now let `t>=1`. Primitivity implies that `u,v` are not both even. If `u`
is odd, the relevant even pair is `(a_4,a_6)` when `v` is odd and
`(a_3,a_5)` when `v` is even. In either case the pair differs by `2u`, of
valuation one, and `a_2=2u` also has valuation one. At `t=1`, at most two
forms can satisfy `s_j=t`, and the other member of the pair contributes at
least one to `D`; hence `Delta_2<=0`. At `t>=2`, `a_2` is below level `t`
and at most one member of the pair can reach level `t`, so `E<=1` and
`Delta_2<=1`. Finally, if `u` is even and `v` is odd, the four forms
`a_3,a_4,a_5,a_6` are odd and only `a_1,a_2` can meet level `t`; their
valuations are consecutive, `s_2=s_1+1`, so `E<=1`.

Therefore

```text
v_(1+i)(G)-v_(1+i)(P) <= 1,
|G|/|P| <= sqrt(2)*5 = 5 sqrt(2).                    (3)
```

The companion checker
[check_sixfactor_ramified_content_sharpening.py](check_sixfactor_ramified_content_sharpening.py)
checks the valuation split exhaustively on a finite representative box and
checks `Norm(G)<=50 P^2` on bounded actual tuples. These computations are
validation only; inequality (3) follows from the case proof above and the
odd-prime argument cited in (1).

## Consequence at `C=1/2`

The exact chord/content identity from
[six_factor_primitive_content_bound.md](six_factor_primitive_content_bound.md)
gives, after primitive division,

```text
normalized arc span >= 2 sqrt(|P|/|G|).
```

Combining this with (3), every arc containing five distinct points in the
family satisfies

```text
normalized arc span >= 2/sqrt(5 sqrt(2))
                       = 0.752120... > 1/2.         (4)
```

Thus no moving choice of integer `u,v,T` in this six-factor
family can produce a five-point primitive tuple on an arc with normalized
length at most `1/2`. In particular there is no nontrivial outer/interior
threshold `K` to grow in that range.

If `T=0` and `P!=0`, every row is `P` or `-P`, so there are at most two
distinct points. If `P=0`, the polynomial chord identity from
[the original family note](five_point_affine_shape_cubic_family.md) remains
valid:

```text
Norm(w_1-w_4) Norm(w_2-w_4)(T^2+4u^2)
 =64 P^2 Norm(w_0)(T^2+u^2).
```

When `T^2+4u^2>0`, its zero right side forces a row collision. When
`T=u=0`, all five products vanish. Thus neither boundary case permits
five distinct points.

For rational parameters, clear a common denominator and divide by
`gcd(u,v,T)`. The resulting primitive projective triple gives the same
primitive Gaussian tuple up to a unit, so (4) covers all rational
clearings as well.

Root and Luna independently derived the ramified case split; Astra
independently audited all cases, the odd-prime combination, and the two
boundary cases. Root replayed the permanent checker successfully on
1,820,078 valuation cases and 174,720 signed tuples. These finite checks
support the displayed all-parameter proof, not a claim about arbitrary
endpoint configurations.
