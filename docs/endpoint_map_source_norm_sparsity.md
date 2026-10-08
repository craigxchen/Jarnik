# Endpoint maps force the source radius into a sparse arithmetic set

The [small-parameter reconstruction](endpoint_map_small_parameter_reconstruction.md)
has a consequence for the source radius alone. The radii that can support
a nonconformal map between two large endpoint configurations form a sparse
set of integers. This gives infinitely many actual primitive circle norms
where such maps are impossible, without finding or estimating a target
radius.

The conclusion concerns maps with many retained points. It does not prove
that an arbitrary endpoint source has a map, or that the excluded circles
have many endpoint points. In particular it does not improve the current
unconditional point-count bound.

## 1. Norm sparsity, uniform in the retained point count

For an integer `M>=128`, let `A_M(X)` be the set of integers `1<=N<=X`
for which there exists a primitive `m`-point source with `m>=M`, least
squared radius `N`, and a rational nonconformal half-angle projective map
to a full primitive target of least squared radius `N'`, such that

```text
Theta N^(1/4) <= 2,       Theta' (N')^(1/4) <= 2.
```

Every label is retained. No ordering of `N,N'` is imposed. Put

```text
beta_M = 1/4 + 8/M + 23/(M-4).
```

Then, for every `epsilon>0`,

```text
#A_M(X) <= C_epsilon X^(beta_M+epsilon),              (1)
```

with a constant independent of `M` and `X`. For example,

```text
beta_128 = 247/496 < 1/2.
```

As `M` grows, the exponent approaches `1/4`. This is an arithmetic
restriction on the possible source norms, rather than a further count
of maps for a fixed source.

## 2. The source divides a controlled quartic value

Retain the notation `E,q,h,u,D_m,rho` from the reconstruction theorem.
Let

```text
K(a,b,d) = (b^2+(a+d)^2)(b^2+(a-d)^2).
```

At the actual anchor selected by that theorem, the map has primitive
triangular representative `((a,b),(0,d))`, with `a,d>0` and `K>0`.
The clipped conductor `Q=N/E` divides its discriminant `K`. Thus

```text
N | E K(a,b,d).                                      (2)
```

This uses the exact clipped conductor of the selected representative.
It is not inferred from the necessary scalar square equation alone.
Nonconformality is essential: when `a=d,b=0`, the discriminant vanishes
and (2) would impose no restriction.

The radius-unordered bounds and the elementary estimates

```text
q<=8/m,       h<=7/m,       u<=1/2+2/m,       D_m<=2^14
```

give, for `m>=128`,

```text
1 <= E <= 2^8 N^(8/m),
1 <= a,d <= 2^15 N^(7/(m-4)),
|b| <= 2^17 N^(1/4+9/(m-4)).                         (3)
```

The four estimates hold for every `m>=128`, not just the finite range
checked below. Indeed, `k>=(m-4)/4`, `m-2k-1>=(m-2)/2`, and
`F>=m(m-2)/4` give

```text
q <= 4m/((m-4)(m-2)) <= 8/m,
g <= 1/(m-2),       ell <= 9/(4(m-1)),
h <= 9/(4(m-1))+2/(m-2) <= 7/m,
log_2 D_m <= 2+4(m-1)/(m-2)+9m/(4(m-1))+m/(m-2) < 14.
```

For `u`, its exact formula with `L=2floor((m+1)/2)>=m` has a
correction of either `1/L` or `L/(L^2-4)`, both at most `2/m`.

For the shear estimate, use `b^2<=T` and the trace estimate from the
reconstruction proof:

```text
|b| <= sqrt(5) 2^(u/2) D_m (16N)^((u/2+h)/rho).
```

Its exponent of `N` is at most

```text
(1/4+8/m)/(1-4/m) = 1/4+9/(m-4),
```

and its absolute factor is less than `2^17` for `m>=128`.
Similarly `D_m(16N)^(h/rho)<=2^15 N^(7/(m-4))`.
These estimates continue to hold when the original target is the larger
circle: apply the reduction to that larger source and take its triangular
adjugate, as in the reconstruction theorem.

Since all exponents in (3) decrease with `m`, the following explicitly
defined arithmetic superset contains `A_M(X)`:

```text
D_M(X) = union of the positive divisors of E K(a,b,d) that are <=X,
```

where

```text
1<=E<=2^8 X^(8/M),
1<=a,d<=2^15 X^(7/(M-4)),
|b|<=2^17 X^(1/4+9/(M-4)),
K(a,b,d)>0.                                         (4)
```

No circle configuration or target data occurs in (4). One can use this
as a sufficient test for nonexistence: if `N` is not in the analogous
divisor set obtained by putting `X=N`, no nonconformal endpoint map with
at least `M` retained points exists at that source norm.

## 3. Counting the possible source norms

The number of parameter quadruples in (4) is at most an absolute constant
times

```text
X^(8/M + 14/(M-4) + 1/4+9/(M-4)) = X^beta_M.        (5)
```

Also `K<=T^2`, where `T=a^2+b^2+d^2`. The shear exponent in (4) is
at least the diagonal exponent. Consequently, for an absolute constant
`C`,

```text
E K <= C X^(1+8/M+36/(M-4)) <= C X^2.                (6)
```

The last inequality holds uniformly for `M>=128`.

For completeness, the elementary divisor bound says that for every
`eta>0` there is `c_eta` with `tau(n)<=c_eta n^eta`.
Indeed for sufficiently large primes, `e+1<=2^e<=p^(eta e)`; for
each of the finitely many smaller primes, the supremum of
`(e+1)/p^(eta e)` over nonnegative integers `e` is finite.
Multiplying these bounds proves the claim. Apply it to (6) with a
sufficiently small `eta` to bound each divisor count by
`O_epsilon(X^epsilon)`. Equations (2), (5), and (6) prove (1).

For a **fixed exact** number of rows `m`, keeping the exact reconstruction
constants improves the exponent in (1) to

```text
beta_m_exact = q+(u/2+3h)/rho
             = 1/4+61/(4m)+O(1/m^2).                 (7)
```

For example `beta_68_exact=96152107/195629280<1/2`.
The uniform theorem uses the coarser threshold 128 to avoid any need to
combine point-count-dependent constants or verify monotonicity of the
piecewise exact constants.

## 4. Infinitely many excluded norms of actual primitive endpoint sources

Consider the explicit family

```text
N_n = 4n^2+1,       n>=1.
```

These are legitimate least primitive squared circle radii, not arbitrary
integer test values. The Gaussian points

```text
z_minus=2n-i,       z_plus=2n+i
```

have common norm `N_n` and Gaussian gcd one: their common divisor would
divide `2i`, but `N_n` is odd. Hence this two-point realization is already
primitive. Its phase width and endpoint constant satisfy

```text
Theta_n = 2 atan(1/(2n)),
Theta_n N_n^(1/4) <= (4n^2+1)^(1/4)/n <= 5^(1/4) < 2.
```

Thus every `N_n` comes with an explicit actual endpoint source of least
radius `N_n`. The displayed source has only two points.

Since `n -> N_n` is injective and `N_n<=4Y^2+1` for `n<=Y`, theorem (1)
with `M=128` proves, for every `eta>0`,

```text
#{n<=Y : N_n supports a full endpoint map with m>=128}
  = O_eta(Y^(247/248+eta)).                           (8)
```

Choosing `eta<1/248` makes this `o(Y)`. Therefore a density-one subset
of these explicit actual source norms permits no such map, at any source
arc and for any number of retained points at least 128.

There is also a wholly arithmetic, effective excluded class: take the
integers `N_n` failing the divisibility test (4) with `X=N_n,M=128`.
The proof of (1) counts this larger arithmetic exceptional set too, so
the failing class is infinite and has density one within the sequence
`N_n`. This avoids defining the excluded class through the unknown
endpoint-map existence property.

The limitation is deliberate. Equation (8) does **not** assert that any
of these circles contains 128 points in an endpoint arc. If such a source
does occur at a norm failing the divisibility sieve (4), it has no second
full nonconformal endpoint realization. Constructing those large sources,
or proving that all sufficiently large endpoint sources instead lie in
the sparse permitted set, remains unproved. The sparse permitted set can
still contain every source relevant to a hypothetical failure of
uniformity.

The [checker](check_endpoint_map_source_norm_sparsity.py) verifies (2) on
actual mapped tuples, the exact exponent in (7), and all uniform bound
exponents and constants through `m=10000`. Its finite checks supplement
the unrestricted counting proof.
