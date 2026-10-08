# Sparse exceptional directions with an approximate slack background

This is a conditional arbitrary-weight theorem for every nonlinear
repeated-Walsh assignment. It does **not** prove that the condition below
holds for general prime weights. The uniform endpoint count remains open.

Use the notation of `walsh_pair_slack_fixed_family_obstruction.md`:
`M=2^t`, `b>=5`, one distinct physical flip at every row, physical
prime-log total `W0`, common content `D`, and

```
kappa = 4 log C - log 4 - D,
q_a = W_a - 4 F_a/M,
u_d = kappa - qhat(d) >= 0       (d != 0).
```

Fix `0<delta<=1`. Suppose there are a number `c>=0` and an exceptional
set `E` of nonzero row directions such that

```
|E| <= M^(1-delta),
|u_d-c| <= (1/8) log M          for every d outside E.       (1)
```

The approximation is essential: exact constancy at two distinct
directions is incompatible with the actual distinct prime logs.
No restriction is imposed on the values on `E`.

For all sufficiently large `M`, depending on `C,delta`, actual source
phases on an arc of length at most `C sqrt(R)` imply

```
W0 >= (1/32) M h log M,
s = floor(delta t/2),
h = 2^floor(log_2[s/(8(log_2(s+2))^2)]).                    (2)
```

In particular `h >= s/(16(log_2(s+2))^2)`, so this gives the growing
factor `h` for arbitrary nonlinear assignments under (1). The source
radius satisfies `log R^2=W0+D>=W0`.

## Selecting a subspace and a coset using the actual weights

Let `N=2^s`. A uniformly random `s`-dimensional linear subspace `H`
contains each fixed nonzero direction with probability `(N-1)/(M-1)`.
Thus

```
E |H intersect E| <= M^(1-delta)(N-1)/(M-1) = o(1).
```

There is an `H` disjoint from `E`. This avoids the exceptional points,
not their linear span; their rank can be as large as `t`.

Call a row low if its assigned prime-log weight is less than
`(1/2)log M`. The assigned rational primes are all distinct, so there
are at most `sqrt(M)` low rows, without needing a prime number theorem.
There are at most `b(M/N-1)` rows whose global assigned label vanishes
on `H`. Averaging the union of these two bad row sets over the `M/N`
cosets of `H` gives one coset `T` containing at most

```
b(1-N/M)+N/sqrt(M) <= b+1                                (3)
```

bad rows; the last inequality uses `N<=sqrt(M)`.

If `b>h`, distinctness of all `b(M-1)` physical primes already gives
`W0 >= (bM/4)log M`, which is stronger than (2) for large `M`.
Hence assume `b<=h`.

Identify `T` with `F_2^s`. Replace the restricted label at every bad
row by any nonzero functional on `H`. Apply the arbitrary nonlinear
aligned-flat theorem to this patched field, with precisely the `h`
in (2). Its family retains at least `N^(3/4)` rows. Discard every
retained coset meeting a patched row. At least

```
N^(3/4)-(b+1)h >= N^(3/4)-h(h+1) > 0                    (4)
```

rows remain for large `M`. Every surviving flat `P=x0+V` is an actual
global aligned flat, with `V<=H`, a nonzero restriction `ell`, and

```
Fp >= (h/2) log M.                                      (5)
```

Patching low rows as well as zero restrictions makes (5) pointwise:
no surviving flat contains a low assigned prime.

## The constant background costs only average slack

Put `aC=max(0,4 log C-log 4)` and `U=sum_(d!=0) u_d`.
The exact pair inversion gives

```
U = W0-4F/M+(M-1)kappa <= W0+(M-1)aC.
```

For large `M`, at least `M/2` directions are outside `E`; (1) implies

```
c <= 2W0/M+2aC+(1/8)log M.                              (6)
```

For the surviving flat define

```
Psi = sum_(d in V minus {0}) (-1)^(ell(d)+1) u_d.
```

The signs in this sum add to **one**, because `ell` is nonzero.
Consequently its constant background contributes `+c`, and not `-c`
or zero. Since `V<=H` avoids `E`, (1) and (6) give

```
Psi <= c+(h-1)(1/8)log M
     <= 2W0/M+2aC+(h/8)log M.                            (7)
```

Use the exact aligned-flat necessary inequality from the pair-slack
note:

```
(M-4)[2Fp-4log(h/2)] <= M Psi+4W0-4kappa-4F.             (8)
```

Also `|kappa|<=W0/(M-1)+aC`. Substituting (5) and (7), and dropping
`-4F`, yields for `M>=5`

```
(M-4)[h log M-4log(h/2)]
 <= 7W0+(2M+4)aC+(Mh/8)log M.                            (9)
```

Since `h` grows and `h=O(t)`, the left side minus the last term on
the right is `(7/8-o(1))Mh log M`. Equation (9) proves (2), with
ample margin, for sufficiently large `M` depending on `C,delta`.

## What this leaves open

Condition (1) is a sufficient condition on the actual pair slack, not
a consequence of its nonnegativity or bounded mean. In particular,
a Ramsey argument would need quantitative control of both the number
of colors and the dimension it retains. Coloring by the first nonzero
coordinate uses `t` colors and avoids a monochromatic plane; this does
not refute a Ramsey theorem with a fixed number of colors. A general
argument still needs to obtain useful control of signed character sums
such as `Psi`, or exploit the simultaneous certificates beyond this
approximate-background case.

The later [Ramsey growth theorem](walsh_arbitrary_weight_ramsey_growth.md)
uses all four-row phase constraints to bound every slack value, and
thereby gives an arbitrary-weight growth improvement without (1).
The stronger rate in (2) here still uses its stated sparse-background
hypothesis.
