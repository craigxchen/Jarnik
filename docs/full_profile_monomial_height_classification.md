# Which multiplicative words fit the endpoint row height?

For a full Boolean profile, the only nonzero orientation-exponent
patterns that fit the original row's modulus, after removing arbitrary
ordinary integer norm factors, are single rows and quotients of two rows,
up to conjugation. This statement allows subpower Gaussian corrections. It
does not cover additive cancellation or division by large Gaussian
factors with uncontrolled phase.

This identifies precisely why multiplying the available thin pair
products does not automatically create the common signed support
required by [the four-witness gap](four_signed_witness_matching_norm_gap.md).
It is a limitation of that construction, not an endpoint count bound.

## 1. An exact integer exponent inequality

Fix `m>=2`, put `r=2^(m-1)`, and for `epsilon in Z^m` define

```text
F(epsilon)=sum_(T subset [m]) |sum_(i in T) epsilon_i|.
```

The empty subset contributes zero. Then

```text
F(epsilon) >= r max_i |epsilon_i|.                 (1)
```

If `F(epsilon)<=r`, the vector is exactly one of

```text
0,       +e_i,       -e_i,       e_i-e_j (i!=j).     (2)
```

Every nonzero vector in (2) has `F=r`. Every other integer vector
has the stronger lower bound

```text
F(epsilon) >= 3r/2.                               (3)
```

For (1), fix `i` and pair the subsets `U` and `U union {i}`,
where `i notin U`. Each pair contributes
`|x|+|x+epsilon_i|>=|epsilon_i|`. There are `r` such pairs.
For equality with a nonzero vector, change the overall sign so
`epsilon_i=1`. Equality forces every other subset sum to be in
`[-1,0]`. Singletons force the other entries to be zero or minus one;
two minus ones violate the condition. This proves (2).

To prove (3), an entry of magnitude at least two gives `F>=2r`.
Otherwise all entries are zero or signed units. A vector supported
on at most two coordinates and absent from (2) has two equal signs,
so `F=2r`. If at least three entries are nonzero, condition a uniform
random subset on all other coordinates. The contribution from these
three coordinates has the distribution `B-c`, where `B` is
`Binomial(3,1/2)` and `c` is an integer. For every integer shift `h`,
`E|B-h|>=3/4`; the minimum is at `h=1,2`. Thus
`F/2^m>=3/4`, proving (3).

## 2. The Gaussian norm budget, including corrections

Take a full extracted profile

```text
P_i=K_i product_(T containing i) H_T,
log N(H_T)=w+o(w),        log N(K_i)=o(w),
```

with `m` fixed, `w` tending to infinity, and pairwise rational-prime-
disjoint conjugate-primitive core blocks. For a multiplicative word
in rows and conjugate rows, let `epsilon_i` be its number of `P_i`
factors minus its number of `bar(P_i)` factors. Balanced pairs
`P_i bar(P_i)=N(P_i)` are ordinary integers and do not affect its
Gaussian orientation imbalance.

At a split prime in `H_T`, that imbalance is the core depth times
`sum_(i in T) epsilon_i`, plus the correction contribution.
Dividing by an ordinary rational number changes the valuations at
the two conjugate primes equally. Any nonzero Gaussian-integer
output must have total valuation at least the absolute imbalance.
Allow also one further Gaussian rational multiplier whose numerator
and denominator have total log norm `o(w)`. Summing the primewise
triangle inequalities gives

```text
log N(output)
 >= sum_T |sum_(i in T) epsilon_i| log N(H_T)
    -sum_i |epsilon_i| log N(K_i)-o(w).             (4)
```

This bound retains corrections supported on core primes; it does not
delete an entire prime-power block because of a small overlap.
Prime factors outside the cores can only increase the output norm.

Each original row has norm `exp(rw+o(w))`. If the output also has
log norm at most `rw+o(w)`, (1) and (4) first bound
`max|epsilon_i|`: writing the uniform core error as `eta(w)w`
and the maximum correction log norm as `kappa(w)w`, the right side
of (4) is at least

```text
max|epsilon_i| [r(1-eta(w))-m kappa(w)]w-o(w),
eta(w),kappa(w) -> 0.
```

Consequently `max|epsilon_i|<=1` for all sufficiently large `w`.
There are now finitely many exponent vectors. Equations (3)--(4)
exclude every vector outside (2), since its core log norm is at
least `3rw/2-o(w)`.

For a signed unit vector the surviving core is one original row
or its conjugate. For `epsilon=e_i-e_j`, rational norm cancellation
gives the familiar pair product

```text
Z_ij=P_i bar(P_j)/G_ij,
G_ij=product_(T containing i,j) N(H_T).
```

Its support is precisely the subsets containing exactly one of
`i,j`. These supports differ for different unordered pairs, and from
the original row supports. The zero exponent case has only subpower
orientation data. It may leave a large ordinary real factor, but it
cannot produce a full signed core of log norm `rw+o(w)` with a
subpower correcting multiplier: such a core has leading orientation
imbalance, while the output does not.

Thus multiplicative words at the original endpoint modulus do not
enlarge the available support classes. A large Gaussian division can
alter this conclusion algebraically, but it also rotates the output;
its imaginary height requires a separate argument. Additive
combinations likewise lie outside (4), since their valuations may
increase by cancellation.

The [exact checker](check_four_signed_witness_parity_gap.py) enumerates
bounded integer exponent vectors for `m=2,...,5` as a supplemental
check of (1)--(3). The proof above covers every integer vector and
retains the full asymptotic correction budget.
