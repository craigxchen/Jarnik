# A conditional divisor bound for critical triangular maps

This note audits the divisor consequence of the [reciprocal-stretch
grid](mobius_reciprocal_stretch_grid.md) against the [finite-sample critical
reduction](mobius_finite_sample_critical_form.md). It bounds the number of
points only when a nonconformal endpoint-to-endpoint map satisfying the
critical reduction exists. It is not a radius-independent point-count
theorem for arbitrary endpoint configurations.

The subsequent [squareclass argument](reciprocal_squareclass_runge_growth.md)
improves this conditional growth rate to
`O(sqrt(log R/loglog R))`. The explicit conductor bounds proved below
remain inputs to that stronger result.

Let a primitive source tuple of `m` points have least squared radius `N`,
and let its image have least squared radius `N'<=N`. Assume both endpoint
constants are at most two and `m>=8`. Use the actual anchor selected by the
critical reduction and its primitive triangular integer representative
`M=((a,b),(0,d))`. Conjugate the whole target if needed so `d>0`, and put
`delta=ad>0`, `h_m=g_m+b_m`, and `Q=Nclip`. The reduction supplies

```text
delta <= D_m^2 N^(2h_m),
Q >= 2^(-m q_m) N^(1-q_m),
h_m = 13/(4m)+O(m^(-2)),       q_m = 4/m+O(m^(-2)).
```

The clipped-conductor theorem gives `Q|N` and `Q|N'`. The projection
content proof gives `Q|g` and `Q|g'`, where `g,g'` are the ordinary
coordinate gcds of the source and inverse target projection normals.
Consequently every normalized squared stretch `x_i>0` produces integers

```text
n_i=2 delta (N/Q) x_i,
n_i'=2 delta (N'/Q) /x_i,
n_i n_i'=H=4 delta^2 N N'/Q^2.
```

The factor `2` is required by the physical projection identity. Both
`n_i,n_i'` are positive integers because `Q` divides the two projection
normals coordinatewise. Each `n_i` occurs at most twice, since a nonzero
linear projection of a circle has at most two points at one value. Hence

```text
H is a positive integer,                 m <= 2 tau(H).       (1)
```

If `Q=N`, then `Q|N'<=N` forces `N'=N`, so `H=4delta^2`. In particular
`delta=1` gives `m<=2tau(4)=6`.

The critical bounds imply

```text
H <= 4 D_m^4 2^(2m q_m) N^(4h_m+2q_m),                 (2)
4h_m+2q_m = 21/m+O(m^(-2)).
```

There is a uniform coarse version. The exact estimates recorded in the
parabolic-anchor audit give, for `m>=256`, `h_m<=7/m`, `q_m<=8/m`, and
`D_m<=2^14`. Thus

```text
H <= 2^74 N^(44/m).                                     (3)
```

The following elementary divisor estimate is enough to invert (1)--(3):

```text
log tau(n) <= (log 2+o(1)) log n/loglog n               (n -> infinity).
```

For a proof without a prime-counting theorem, write `y=log n` and split
the prime factors of `n` at `R=y/(log y)^3`. At primes `p<=R`, there are
at most `R` possibilities, and each exponent is at most `y/log 2`, so
their contribution to `log tau(n)` is at most
`R log(1+y/log 2)=O(y/(log y)^2)`. At primes `p>R`, use
`log(a+1)<=a log 2` and `sum a log p<=y`. Their contribution is at most
`(log 2)y/log R=(log 2+o(1))y/log y`.

Write `W=log N`. If `m>=100 W/(log W loglog W)` along arbitrarily large
`W`, then (3) gives

```text
log H <= 74 log 2 + (44/100) log W loglog W.
```

For unbounded `H`, the divisor estimate yields
`log(2 tau(H)) <= ((44/100)log 2+o(1))log W`, while the assumed lower
bound on `m` gives `log m>=(1+o(1))log W`. Since
`(44/100)log 2<1`, this contradicts (1). If `H` stays bounded, (1)
contradicts the same lower bound directly. The finitely many `m<256`
are harmless. Therefore every map in this stated admissible class obeys

```text
m = O(log N/(loglog N logloglog N))                    (N -> infinity).
```

This result uses the map, the two endpoint assumptions, `N'<=N`, and the
selected actual-anchor triangular reduction. It does not establish that
such a map exists for every endpoint configuration, and it does not make
the determinant or conductor loss uniformly bounded.
