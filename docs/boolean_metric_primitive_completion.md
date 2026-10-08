# Primitive Gaussian completion with only residual-height loss

The joint norm, Pluecker, and metric equations can supply a Gaussian
integer realization by the
[metric realization criterion](boolean_norm_metric_realization_criterion.md).
That realization need not initially be conjugate-primitive or exhibit
the prescribed core factors in one orientation. This note proves that
both properties can be recovered after only subpower changes to the
core norms and correcting factors. The needed common rotation preserves
angular differences; it need not preserve small imaginary coordinates
relative to the integer real axis.

There is no strict height gain. The remaining fixed-metric rotation
problem is described in the
[common-rotation classification](gaussian_common_rotation_classification.md).

## Statement

Fix `m>=2`. Suppose Gaussian integers `P_i` have positive odd norms
`N_i` and nonzero determinants, with exact factorizations

```text
N_i=k_i product_(T containing i) n_T,
Delta_ij=Im(bar(P_i) P_j)
        =t_ij product_(T containing i,j) n_T,
k_i in Z_(>0),           t_ij in Z\{0}.
```

Here the `n_T` are positive odd integers with pairwise disjoint
rational-prime support, all supported on split primes. Unit core
factors are allowed. Set

```text
K=product_i k_i,              B=product_(i<j)|t_ij|.
```

There exists a common rational Gaussian rotation `u`, with
`Norm(u)=1` and all `uP_i` Gaussian integral, such that, on putting

```text
g_i=gcd(|Re(uP_i)|,|Im(uP_i)|),
R_i=uP_i/g_i,
```

the rows `R_i` are conjugate-primitive and

```text
product_i g_i divides K B.                                (1)
```

There also exist oriented Gaussian blocks `H'_T`, with
`n'_T=Norm(H'_T)` dividing `n_T`, for which

```text
R_i=K'_i product_(T containing i)H'_T,
Im(bar(R_i) R_j)=t'_ij product_(T containing i,j)n'_T,
K'_i in Z[i],             t'_ij in Z\{0},
L=sum_T log(n_T/n'_T) <= 2 log B+log K,                    (2)
log Norm(K'_i) <= log k_i+L,
log|t'_ij| <= log|t_ij|+L.                                (3)
```

Every `K'_i` is itself conjugate-primitive, and no orientation opposite
to an incident retained block can occur in it. The retained blocks
still have disjoint rational-prime support.

Consequently, for fixed `m`, if every original core has
`log n_T=w+o(w)` and `log K+log B=o(w)`, then all retained cores
have `log n'_T=w+o(w)`, while all the new correcting factors, row
contents, and determinant residuals have logarithmic height `o(w)`.

## 1. The local common-core lemma

Fix a split prime `p=pi bar(pi)` in the unique core `n_T`, and
write `e=v_p(n_T)>0`, `kappa_i=v_p(k_i)`. For a Gaussian integer
row write

```text
a_i=v_pi(P_i),          b_i=v_bar(pi)(P_i).
```

For `i in T` we have `a_i+b_i=e+kappa_i`. If `i,j in T`, then
`v_p(Delta_ij)>=e`. The two Gaussian valuations of
`bar(P_i)P_j` are

```text
alpha=b_i+a_j,          beta=a_i+b_j.
```

Both `alpha,beta` are at least `e`. Indeed, if they differ, the
valuation of the difference between this product and its conjugate is
their minimum. Division by `2i` costs nothing because `p` is odd, so
that minimum is `v_p(Delta_ij)>=e`. If they are equal, each is
`e+(kappa_i+kappa_j)/2>=e`.

Set

```text
A_T=min_(i in T)a_i,       B_T=min_(i in T)b_i.
```

Choose rows attaining these two minima. If they are different, the
preceding product inequality proves `A_T+B_T>=e`; if they coincide,
their row norm proves the same conclusion. The upper bound follows
by comparison with every row norm. Thus

```text
e <= A_T+B_T <= e+min_(i in T)kappa_i.                    (4)
```

This also proves (4) for a singleton `T`. It does not assert that the
Gaussian common factor is already oriented: both `A_T` and `B_T`
may be positive.

## 2. Rotate only the full-core part

If `T=[m]`, the common Gaussian gcd of all the rows has local exponents
`A=A_T`, `C=B_T`. Multiplying every row by
`(pi/bar(pi))^C` replaces those exponents by `A+C,0` and keeps
every row integral. Perform this operation at each prime in the full
core, and no other primes. The resulting common multiplier is the
claimed rational rotation `u`.

Write `A+C=e+delta`. Equation (4) gives
`0<=delta<=min_i kappa_i`. After removing that common factor, the
local norm of row `i` has exponent `kappa_i-delta`. Therefore the
ordinary rational content exponent `c_i` of the rotated row satisfies

```text
c_i<=kappa_i,                 p in the full core.        (5)
```

After dividing each row by its ordinary content, a common factor in
the selected orientation still has depth at least

```text
e'_full=max(0,e-sum_i kappa_i).                           (6)
```

Rotations at these full-core primes do not change any valuation at a
different rational prime.

## 3. Proper-core orientation and row contents

Suppose `T` is proper. Define

```text
rho_ij=v_p(t_ij),
R_p=sum_(i in T,j notin T)rho_ij.
```

Each crossing pair occurs once in this sum. Since `p` belongs to
no core shared by a crossing pair, its determinant valuation is exactly
`rho_ij`. If `c_i=min(a_i,b_i)` is the ordinary row-content
exponent, the rational integer `p^c_i` divides `P_i` and hence every
determinant involving that row. Consequently

```text
c_i<=min_(j notin T)rho_ij   (i in T),
sum_(i in T)c_i<=R_p.                                    (7)
```

For an outside row, `a_i+b_i=kappa_i`, so
`c_i<=kappa_i/2`.

The rational content of the `T`-common Gaussian factor is
`min(A_T,B_T)`. It is at most every `c_i` with `i in T`, and
hence at most `R_p`. By (4), its larger oriented exponent is at
least `e-R_p`. Dividing the incident rows by their respective
ordinary contents loses at most a further `R_p`. Therefore one
orientation of `p` divides all the primitive incident rows to depth

```text
e'_T=max(0,e-2R_p).                                      (8)
```

The orientation can depend on `p`; each Gaussian block is formed by
collecting these compatible oriented prime powers. No requirement of
whole-block conjugation is imposed.

At a prime outside every core, ordinary row content is bounded by
half the valuation of its norm `k_i`. Combining this fact, (5), and
(7), prime by prime, proves (1). The prime `2` is absent from all
row norms by hypothesis. Dividing by the ordinary contents therefore
leaves odd-norm Gaussian integers with coprime integer coordinates,
which are exactly conjugate-primitive Gaussian integers.

## 4. The total cost and new residuals

Use the exponents (6) and (8) to define `H'_T`. They divide every
incident primitive row, with one common orientation at every prime.
Their norm losses obey

```text
sum_(proper cores T) log(n_T/n'_T)
 <=2 sum_(proper-core primes p)R_p log p <=2 log B,
log(n_[m]/n'_[m]) <=log K.
```

The core supports are disjoint, so no prime is charged to two core
positions. This proves (2).

Gaussian division now defines integral `K'_i`. A common Gaussian
factor divides two rows only if its norm divides their determinant,
which proves integrality of the displayed `t'_ij`. Nonzero determinants
are preserved. More explicitly,

```text
Norm(K'_i)=k_i product_(T containing i)(n_T/n'_T)/g_i^2,
t'_ij=t_ij product_(T containing i,j)(n_T/n'_T)/(g_i g_j).
```

Taking logarithms and dropping the nonpositive content terms proves
(3). Conjugate-primitivity of `R_i` immediately gives the corresponding
assertions about `K'_i` and its incident orientations.

## Scope

The metric realization theorem followed by this result supplies an
actual conjugate-primitive Boolean Gaussian factor profile at the same
leading height. It does not supply the endpoint's small nonzero
imaginary coordinates. The transformation has the exact direction law

```text
R_i/bar(R_i)=(u/bar(u)) P_i/bar(P_i),
```

so angular clustering is preserved, and ordinary row-content removal
does not rotate individual rows. But a common rational rotation can
move the cluster away from the integer real axis. The choice of a
norm representation of the common Gaussian gcd that restores small
imaginary coordinates remains an additional simultaneous Diophantine
condition. Neither this primitive completion nor the metric equations
alone establish it.

The [exact valuation checker](check_boolean_metric_primitive_completion.py)
exhausts small split-prime valuation patterns and checks additional
larger patterns. It verifies the local common-core bound, row-content
bound, retained oriented depths, and loss budgets. The proof is
primewise and applies to arbitrary valuations; the computation is a
supplement, not the source of that generality.
