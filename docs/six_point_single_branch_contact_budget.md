# The single-branch contact budget

Let `X/Q` be the geometrically integral intersection of ten quadrics in
`P^11`, so `deg(X)=1024`. Let `B` be a branch point of residue degree `D`.
A homogeneous auxiliary form followed only by coefficient-field norm and
a nonzero-integer lower bound cannot overcome the full branch orbit when
`D=12288`.

Indeed, let `F` be a degree-`d` homogeneous form over a number field `E` of
degree `e`, and assume `F|_X` is nonzero. Here `d` and `e` are fixed
independently of the endpoint height, and the coefficient heights are
`U^O_(d,e)(1)`. Its coefficient norm

```text
S = Norm_E/Q(F)|_X = product_(sigma:E->Qbar) F^sigma|_X
```

is a nonzero rational section of `O_X(ed)`. Define its total contact at B by

```text
r = ord_B(S) = sum_sigma ord_B(F^sigma).               (1)
```

Thus `r` includes all conjugate contacts, not only the contact of the
original form. Rationality gives the same order `r` at every one of the `D`
conjugates of `B`. Since the zero divisor of `S` has degree `1024ed`,

```text
D r <= 1024 e d.                                      (2)
```

At an endpoint approaching B, the local branch parameter has the upper
bound

```text
|zeta| <= U^C T^-10.
```

After absorbing the fixed-degree coefficient and chart heights, evaluation
of the norm section gives

```text
|S(z)| <= U^C T^(ed-10r),
ed-10r >= ed(1-10240/D).                              (3)
```

For `D=12288`, the exponent in (3) is at least `ed/6>0`. Clearing
polynomially bounded denominators and using only that a nonzero integral
value has absolute value at least one therefore yields no decaying power of
`T` for any such fixed `d` and `e`. The proof does not assume that `E`
contains `Q(B)`.

The pointed-cover theorem in
[six_point_pointed_cover_four_torsion.md](six_point_pointed_cover_four_torsion.md)
forces `D=12288` for normalized twists with a rational lift when the branch
group is `S_12` and the discriminant squareclass is not `-1`. In those cases
the present budget rules out obtaining a negative target-height exponent
from single-branch homogeneous contact plus the integer norm, at every
fixed auxiliary degree. It is not restricted to the linear osculating
hyperplane.

This is only an exponent budget. If `F|_X=0` or its endpoint value is zero,
there is no nonzero norm to which the lower bound applies. Rational
functions with poles, finite-place divisibility, auxiliary divisors using
other orbits, and methods not reducible to this one norm comparison are
outside the statement. It is neither an endpoint lower bound nor a barrier
to all methods.
