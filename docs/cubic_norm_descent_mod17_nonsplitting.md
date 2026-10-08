# A mod-17 obstruction to splitting the cyclic cubic quartic family

The one-parameter constant-imaginary norm construction over the cyclic
cubic field below never splits into four linear factors over its original
coefficient field. This is a complete result for this specified family,
not an obstruction to other cubic constructions or to the four-row
Boolean profile.

Let `alpha^3-3alpha+1=0`, `K=Q(alpha)`, and put

```text
w=2-alpha^2,       v=alpha^2+2alpha-2,       delta=3,
u=2r/(1+3r^2),     s=(1-3r^2)/(1+3r^2),     h=3u^2,
k=s/(4h),          r in Q\{0}.
```

The quartic from the cubic trace-form construction is

```text
a=-(t^3+3t)/2+k(t^2+1)^2,
b=(t^2+1)(kt^2-t/2+k(1-2h)),
c=(t^2+1)(2kut-(h+1)/(2*3*u)),
U=a+bw+cv,           Y=[k(1+w)]^(-1),          P=Y(U+i).
```

It is monic, has the fixed root `i`, has a nonzero constant imaginary
part, and `Im Norm_(K(i)/Q(i))(P)=-4 Norm_(K/Q)(Y)`. The exact checker
[check_cubic_norm_descent_mod17.py](check_cubic_norm_descent_mod17.py)
independently computes the field determinant for `r=1/2` here and in the
pure cubic field. The theorem proved below concerns the residual monic
cubic `F=P/(t-i)` for **every** nonzero rational `r`:

> **Theorem.** `F` has nonsquare discriminant in `K(i)` or has a
> non-split reduction at a place above 17. In particular, `F` does not
> split completely over `K(i)`.

## The residual cubic

Let `A=2h/s` and `B=2u/s`. Direct division by `t-i` gives

```text
F(t)=t^3+(b0/a0)t^2+(c0/a0)t+d0/a0,
a0=1+w,
b0=i*a0-A*a0+2u*v,
c0=1+w*(1-2h)-i*a0*A+v*B*(-h-1+i*s),
d0=-2A+i*(1+w*(1-2h)-v*B*(h+1)).                 (1)
```

One way to check (1) is to note that
`(U+i)/(t-i)=-1+(t+i)R(t)` with `R` quadratic, and then divide by
its leading coefficient `k(1+w)`. All denominators in (1) are nonzero
in `Q`: `r!=0`, `1+3r^2>0`, and `1-3r^2!=0` for rational `r`.

The polynomial `x^3-3x+1` has simple roots `7` and `13` modulo 17,
and `x^2+1` has the simple root `4`. Hensel's lemma therefore gives
two embeddings `K(i) -> Q_17`, labeled by `alpha=7,13` and `i=4`
modulo 17. At these two places, respectively,

```text
(w,v,1+w)=(4,10,5), (3,6,4) mod 17.                 (2)
```

## Parameters that are 17-adic units

If `v_17(r)=0`, both `1+3r^2` and `1-3r^2` are units: `-1/3=11`
and `1/3=6` are nonsquares modulo 17. Thus (1) is an integral monic
cubic at each of the two places, with stable degree. At `alpha=7`,
`i=4`, its reduction has nonzero discriminant and fewer than three
roots in `F_17` for every `r mod 17` except `2,3,16`. At the second
place the three exceptions reduce to

```text
r=2:  t^3+5t^2+t+7    =(t-15)(t^2+3t+12),
r=3:  t^3+14t^2+2t+16=(t-6)(t^2+3t+3),
r=16: t^3+14t^2+7t+14=(t-12)(t^2+9t+13),             (3)
```

all modulo 17. The quadratic discriminants are `12,14,12`, all
nonsquares. The `r=16` reduction at the first place has repeated
roots, so the second place is essential there. Formula (1) and a loop
over the sixteen residues, without a search over unbounded rational
parameters, give the finite certificate in the checker. A split monic
integral cubic would reduce to three linear factors, contrary to this
certificate.

## Parameters near zero or infinity at 17

Use the first place (`alpha=7`, `i=4`). If `v_17(r)>0`, set `x=r`;
if `v_17(r)<0`, set `x=1/r`. In either case the coefficients of `F`
belong to `Z_17[[x]]`: the denominators in (1), after substituting
`r=x` or `r=1/x` and canceling powers of `x`, have constant terms
among `1`, `3`, `1+w`, and their products, all 17-adic units. They
tend as `x -> 0` to

```text
F_0(t)=(t-i)(t+i)^2.
```

To retain the double-root splitting, put `t=-i+xz`. Exact expansion
of (1) gives the following leading tangent quadratics:

```text
F(-i+xz)/x^2 ->
  -2i*z^2 -4(1+i)v*z/(1+w) -16*delta/(1+w)
                                      when x=r -> 0,
  -2i*z^2 +4(1-i)v*z/[delta*(1+w)]
             +16/[delta*(1+w)]      when x=1/r -> 0.  (4)
```

Their discriminants, before reduction, are

```text
D_0= 32i [v^2-4delta(1+w)]/(1+w)^2,
D_inf=-D_0/delta^2.                                  (5)
```

At the chosen place, (2) gives `D_0=11` and `D_inf=12` modulo 17;
both are nonsquares. More explicitly, the **full monic cubic**
discriminant satisfies

```text
Disc(F)/r^2     = c_0+O(r),     c_0=-4D_0  =7 mod 17,
Disc(F)*r^2     = c_inf+O(1/r), c_inf=-4D_inf=3 mod 17. (6)
```

in `Z_17`. The exact formal power-series computation in the checker
verifies the zero constant and linear terms and the coefficients
`7,3` of `x^2`. These coefficients also follow from (5): the simple
third root tends to `i`, and the leading full discriminant is
`-4D_0 x^2` or `-4D_inf x^2`. Both `7` and `3` are nonsquares modulo
17. Since `r^2` and `r^(-2)` are squares, (6) makes `Disc(F)` a
nonsquare in `Q_17` for every nonunit parameter. The discriminant of
a completely split cubic is a square, proving the theorem.

This fixed-place obstruction closes the splitting question for this
explicit cyclic one-parameter family. It does not rule out a different
choice of cubic trace-form vectors, another field, or a quartic not in
this family.

It also persists under any nonconstant polynomial composition. If
`f in K(i)[t]` and `P(f(t))` split completely over `K(i)`, then every
root `beta` of `P(f(t))` would lie in `K(i)`. For every root `gamma`
of `P`, the nonconstant polynomial `f(t)-gamma` has a root `beta` over
an algebraic closure, so `gamma=f(beta)` would belong to `K(i)`.
That contradicts the theorem. In particular, a rational real
quadratic `f(t)=t+a(t^2+1)` fixes `i` and preserves the constant
imaginary part, but cannot turn this quartic into a split degree-eight
row over the same field.

## Rotated trace-form bases: a unit-parameter result

There is a bounded extension of the residue test. Rotate the two
trace-zero vectors by rational `p,q` with `p^2+3q^2=1`:

```text
w'=p*w+q*v,          v'=-3q*w+p*v.                    (7)
```

Use `w',v'` in the same quartic formulas. Every rational rotation has
17-adically integral `p,q`, by the anisotropy argument in the
[complete classification](cubic_quartic_complete_parameterization.md).
Suppose now that `r` is a 17-adic unit. The checker enumerates the
18 reductions of `(p,q)` on `p^2+3q^2=1` over `F_17`. For fifteen
rotations, all three embeddings have unit leading coefficient
`1+w'`; for each of the sixteen nonzero `r` residues at least one of
their residual cubics has a non-split reduction at `i=4`.

Three rotations have a zero leading coefficient at two embeddings:

```text
(p,q)=(2,13), (5,3), (10,1) mod 17.
```

For each, the three `(w',v')` pairs across the embeddings are exactly
`(2,0),(16,3),(16,14)`, in different orders. At the regular pair
`(2,0)`, the residual cubic reduces to three linear factors only for
`r=2,5,12,15 mod 17`. Those four cases are eliminated at one of the
two bad-leading places, allowing both local roots `i=4,13`:

| `r mod 17` | `(w',v')` | `i mod 17` | Reduction of `(U+i)/(t-i)` |
|---:|:---:|---:|:---|
| 2 | `(16,3)` | 4 | `5t^2+4` |
| 5 | `(16,14)` | 13 | `2t^2+16t+5` |
| 12 | `(16,3)` | 13 | `2t^2+16t+5` |
| 15 | `(16,14)` | 4 | `5t^2+4` |

The two displayed quadratics are irreducible modulo 17: their
discriminants are respectively `5` and `12`, both nonsquares. The
unnormalized cubic has integral coefficients
and is primitive at these places because its quadratic coefficient is
a unit. By Gauss's lemma over `Z_17`, a primitive polynomial that
splits into linear factors over `Q_17` reduces to a product of linear
or constant factors. The irreducible quadratic reductions therefore
exclude splitting despite the leading-degree drop. This proves
nonsplitting for all rationally rotated bases when `r` is a 17-adic
unit. The [local boundary proof](cubic_rotated_mod17_local_escape.md)
exhibits a rational rotation for which every parameter of valuation
at least three splits at all six places above 17. Hence this last
result cannot be extended to all nonunit parameters at the same prime.

This prime-17 unit obstruction does not persist at every split prime.
At `197`, take the rational rotation
`p=-971/973`, `q=36/973` (so `p^2+3q^2=1`) and `r=102`.
The three roots of `alpha^3-3alpha+1` modulo 197 are `34,169,191`;
the two roots of `i^2+1` are `14,183`. All parameter denominators and
all six leading coefficients are units. At the six places, in the
order `(alpha,i)=(34,14),(34,183),(169,14),(169,183),(191,14),
(191,183)`, the residual cubics have three distinct roots respectively

```text
(8,72,120), (70,177,178), (3,126,152),
(158,165,183), (24,73,75), (10,63,127) mod 197.
```

The checker verifies the rational rotation identity and all six
coefficient and root lists exactly. This is a fully split **local
reduction** at 197; it does not make the cubic split over `K(i)`.
