# A residual-discriminant bound for the normalized local kernel index

Fix `n,d>=2`, put `m=nd`, and partition the labels into `T` and `O`.
Let `R` be a characteristic-zero DVR with fraction field `F`, and write

```text
z_i=s a_i  (i in T),       z_j=b_j  (j in O),
```

with `s` any nonzero element of `R`, including a unit. The parameters
`a_i,b_j` lie in `R`.
Assume that the directions `z_1,...,z_m` are pairwise distinct in
`F`. Define the *residual discriminant*

```text
D_T(a,b) = (product_{i<k in T}(a_i-a_k))
           (product_{j<l in O}(b_j-b_l))
           (product_{j in O} b_j).
```

Assume also that `D_T` is nonzero; it need not be a unit. (Distinct
directions alone do not exclude an outside value `b_j=0`.) For
`|I|=2n`, let
`g_I=max(0,|I cap T|-n)` and let `Lbar_I=s^{-g_I} det[V|diag(z)V]_I`.
The normalized determinant is polynomial in `s,a,b` over `Z`. Let
`M_R` be the integral multilinear invariant lattice, `N_R` the span
of all `Lbar_I` times complementary products of `d-2` ordinary
`n`-row determinants, and `K_R` the saturated coherent kernel.

**Theorem.** There are a positive integer `A_(n,d)` and an exponent
`H_(n,d)>=0`, depending only on `n,d`, such that for every `T`, every
characteristic-zero DVR, and every tuple above,

```text
A_(n,d) D_T(a,b)^H_(n,d) K_R  subset  N_R  subset  K_R.       (1)
```

More precisely, if `v_R` is normalized by `v_R(pi)=1`, then

```text
length_R(K_R/N_R)
    <= v_R(A_(n,d)) + H_(n,d) v_R(D_T(a,b)).                   (2)
```

Neither constant depends on `v_R(s)`. The theorem is an existence
bound: the proof gives a finite Gröbner-basis procedure to obtain
`A_(n,d),H_(n,d)`, but does not give useful numerical values. When
`s` is a nonunit and `D_T` is a unit, the sharper bound
`C_(n,d)K_R subset N_R` from the
[linked determinant theorem](linked_determinantal_coalition_saturation.md)
still applies; `C_(3,4)=144`.

## Constant rank away from the residual discriminant

Work first over an algebraically closed characteristic-zero field
`k`. Fix values of `s,a,b` with `D_T(a,b)!=0`. The matrix of
normalized generators has the same rank `h(n,d)` at every such point.

If `s=0`, take the DVR `k[[t]]` and replace `s` by `t` while keeping
the residual values constant. Every factor of `D_T` is a unit there.
The linked determinant theorem says the normalized generator lattice
is the saturated coherent kernel. This kernel is a direct summand of
the free invariant lattice over the DVR, so reduction modulo `t`
preserves its rank. Its generic rank is `h(n,d)` by distinct-configuration
generation. Thus the specialized generator matrix at `s=0` has rank
`h(n,d)`.

If `s!=0`, directions within each of `T` and `O` are pairwise
distinct, since their corresponding factors in `D_T` are nonzero.
A direction can therefore occur at most twice, once in each group.
As `n>=2`, the maximum multiplicity is at most `n`. The
[expected-height calculation](distinct_configuration_determinantal_kernel_generation.md#6-exact-height-when-directions-repeat)
gives height `m-2n+1` for the maximal-minor ideal of
`[V|diag(z)V]`. Eagon--Northcott is exact at that height. Extracting
row multidegree `(1,...,1)` and the `det^d` isotypical part is exact
in characteristic zero. Its terms have dimensions independent of
the values of `z`, so the Euler characteristic of this resolution
shows that the invariant part of the maximal-minor ideal has the
same dimension `h(n,d)` as at distinct directions. It is precisely
the span of the determinant generators with invariant complementary
coefficients. The complementary invariants are spanned by products
of ordinary determinants. As `s` is invertible, normalization does
not change this span. This also proves rank `h(n,d)` when there are
cross-group coincidences.

For `(n,d)=(3,4)`, the
[resolution calculation](coherent_kernel_eagon_northcott_resolution.md)
gives `h(3,4)=341`.

## Universal maximal-minor certificate

The integral multilinear invariant lattice has an integral
determinant-tableau basis, so choose it once and write all normalized
generators as the columns of a finite matrix

```text
G_T(s,a,b) over B_T=Z[s,(a_i)_(i in T),(b_j)_(j in O)].
```

Put `h=h(n,d)`. The constant-rank result says that, over
`Qbar`, the common zero locus of all `h`-by-`h` minors of `G_T` is
contained in `D_T=0`. Hilbert's Nullstellensatz therefore gives

```text
D_T^H_T in I_h(G_T) tensor Q
```

for some `H_T`. Clear rational denominators to obtain a positive
integer `A_T` and the polynomial identity

```text
A_T D_T^H_T = sum_J c_J(s,a,b) det(G_T)_J,
                     c_J in B_T,                              (3)
```

where `J` runs through `h`-by-`h` submatrices. There are only
finitely many coalitions `T`; multiplying the `A_T` and enlarging
the exponents gives common `A_(n,d),H_(n,d)` for all of them.
This is a certificate about the *ideal of all maximal minors*,
not about one selected large minor.

Specialize (3) to `R`. The generator map has rank `h` over `F` by
the distinct-configuration theorem. Its Smith exponents are
`u_1,...,u_h>=0`. Consequently

```text
I_h(G_T(R))=(pi^(u_1+...+u_h)),
length_R(K_R/N_R)=u_1+...+u_h.
```

Identity (3) implies `sum u_i <= v_R(A_(n,d))+H_(n,d)v_R(D_T)`,
proving (2). Since every `u_i` is at most their sum, it also proves
the annihilator inclusion (1). The field-linear generation theorem
identifies the fraction-field span of `N_R` with the coherent kernel,
so its saturation is indeed `K_R`.

## Scope for a global radius argument

If an application controls the sum of local residual-discriminant
valuations, (2) controls the product of these *local presentation
indices* with no dependence on coalition depth `v_R(s)`. It does
not bound the Archimedean covolume of the saturated global kernel.
In particular, a large minor of the redundant generator matrix
cannot be used as a lower bound for that covolume: the local index
is measured by the gcd of **all** maximal minors, as in (3).
