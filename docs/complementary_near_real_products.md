# No forced complementary pair products in the full fair profile

The proposed product lift is exact: if (A_{ab}) is the oriented Gaussian
pair factor, then

```text
Z=(N/P_ab) A_ab^2 = z_a conjugate(z_b),
```

so searching for (A_{ab},A_{cd}) with
`Norm(A_ab)*Norm(A_cd)=N` asks whether two pair supports complement every
prime layer.  The full equal-cut profile supplies a precise obstruction to
forcing such pairs from valuation data alone.

## Direct obstruction for actual pairs when `C<=2`

For an actual distinct pair in the fixed-unit cluster, let `Delta_ab` be its
angular separation and choose the oriented factor `A_ab` with
`arg(A_ab)=Delta_ab/2` modulo `pi`.  Since `A_ab` is a nonreal Gaussian
integer, `|Im(A_ab)|>=1`; hence

```text
1 <= sqrt(P_ab) sin(Delta_ab/2),
P_ab >= csc^2(Delta_ab/2).                            (1)
```

In the short-arc regime `0<Delta_ab<=Delta=C N^(-1/4)<pi`, strictness of
`sin(t)<t` gives

```text
P_ab > 4/Delta_ab^2 >= (4/C^2) sqrt(N) >= sqrt(N).     (2)
```

Therefore any two actual pair norms satisfy
`P_ab P_cd>N`, independently of their valuation supports. Exact
complementary products are already impossible for `C<=2`; no pairing or
support argument is needed. For a larger fixed arc constant, partitioning
into subarcs with constant at most two only gives this obstruction within
each subarc and does not create complementary products across subarcs.

## Full fair construction (ideal valuation model)

Take `M>=3` source labels and one positive layer of equal formal weight
`w` for every nontrivial unordered cut `T|T^c`.  There are

```text
Q=2^(M-1)-1
```

cuts, so `log N=Qw`.  For an unordered source pair `e={a,b}`, let
`chi_e(T)=1` when `T` separates `a` and `b`, and zero otherwise.  Its pair
factor has logarithmic norm

```text
log Norm(A_e) = w sum_T chi_e(T) = 2^(M-2) w
              = (Q+1)w/2.                            (3)
```

The count is independent of `e`: exactly `2^(M-2)` unordered cuts separate
any fixed pair.  Since

```text
(Q+1)/2 <= Q/2 + Q/M       (M <= 2Q),                  (4)
```

every one of the `binom(M,2)` pairs lies in the nearcritical band used by the
Plotkin argument, `log Norm(A_e)<=log N/2+log N/M`.  Thus this profile has
quadratically many nearcritical supports.

Exact positive equal weights on disjoint actual prime supports are
impossible, so this paragraph is an ideal valuation model.  For an actual
near-fair realization with fixed `M`, choose distinct split Gaussian primes
`pi_T` and exponents `s_T` so that

```text
s_T log Norm(pi_T) = w + O_M(1).
```

Orient the prime power `pi_T^(s_T)` on `T` and its conjugate on `T^c`, and set

```text
z_a = product_T pi_T^(s_T 1_{a in T})
             conjugate(pi_T)^(s_T 1_{a not in T}).
```

Every `z_a` then has the same norm `N`, and the quotient `z_a/z_b` uses
precisely the layers with `chi_{ab}(T)=1`.  For fixed `M`, the `O_M(1)` errors
are eventually dominated by the positive nearcritical margin
`(Q/M-1/2)w`, so all pairs remain in the nearcritical window as `w` tends to
infinity.  This is an arithmetic realization of the valuation obstruction;
the Gaussian prime arguments remain uncontrolled.

## Why no two products have norm (N)

For any two source pairs `e,f`, the exponent of the prime layer indexed by
`T` in `Norm(A_e A_f)` is

```text
w(chi_e(T)+chi_f(T)).                                  (5)
```

In the ideal equal-weight model, (3) gives

```text
log Norm(A_e A_f) = (Q+1)w > Qw = log N.                (6)
```

For the actual near-fair prime-power realization above, the same comparison
is `(Q+1)w+O_M(1)` versus `Qw+O_M(1)`, so the product norm is strictly larger
than `N` for sufficiently large `w`.

Distinct source pairs have distinct cut-support vectors `chi_e`; because the
prime supports are disjoint, their ordinary pair norms are consequently
distinct by unique factorization.

There is also a support distinction: if `e` and `f` are distinct, there is
always a cut separating both pairs:
take the common endpoint when they share one, or a cut containing one
endpoint from each pair when they are disjoint.  Hence some `T` has
`chi_e(T)+chi_f(T)=2`, so

```text
Norm(A_e A_f) != N.                                    (7)
```

If `e=f`, every separating cut has coefficient two as well, and equality is
again impossible.  Exact equality `Norm(A_e)Norm(A_f)=N` would require
`chi_e(T)+chi_f(T)=1` for every cut, which the full cut family forbids.

This binary profile shows that the quadratic nearcritical count does not
force complementary pairs or exact norm-`N` products.  The conclusion is
about this disjoint binary valuation profile (and its near-fair
prime-power realizations); arbitrary nested layers can have additional
cross-layer cancellation and are not covered by the exact support argument.

## Geometric scope

The valuation construction is compatible with the pair quotient identities:
for each cut layer, assign one Gaussian orientation to `T` and the conjugate
orientation to `T^c`; then `z_a/z_b=A_{ab}/\overline{A_{ab}}`.  To turn the
profile into an actual shrinking-arc family one would additionally need to
choose the Gaussian prime arguments so that all source phases lie in the
required arc.  That analytic/angular realization is separate.  The exact
conclusion here is that quotient compatibility, equal fair cut weights, and
the nearcritical norm inequality alone cannot supply the complementary
products sought by the proposed route.

If a future geometric argument forces a special pairing of source labels or
imposes phase relations beyond these valuation identities, it would be using
additional information absent from this full-fair obstruction.
