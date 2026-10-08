# Reciprocal fourth moments concentrate on each critical support

Let `m=2s+1>=5`, and let `z_1,...,z_m` be nonzero real numbers with
pairwise distinct absolute values satisfying
`sum_j z_j^(2r+1)=0` for `0<=r<s`. Define

```
p_l=sum_j z_j^(-l),       rho=p_4/p_2^2.
```

Then

```
rho > (m+8)/(9m),
p_2^2/p_4 < 9m/(m+8) < 9.                            (1)
```

Thus the effective support measured by reciprocal squares is uniformly
bounded on every critical cross support. This is a necessary condition
on the full odd moments, not a lattice-circle bound or a solution of the
eight-row problem. It complements the forward moment inequalities in
the [cross-support Hankel note](critical_moment_cross_support_hankel.md).

## Exact identities and a three-vector Gram determinant

The root polynomial is `R(X)=F(X)+r`, where `F` is odd and `r!=0`.
Consequently the elementary symmetric functions of `y_j=1/z_j` satisfy
`e_2(y)=e_4(y)=0`. Newton's identities give

```
p_1^2=p_2,       4p_1 p_3=p_2^2+3p_4.                (2)
```

The vectors `1,y,y^2` are linearly independent: a nonzero quadratic
cannot vanish at the `m>=5` distinct values `y_j`. Their Gram matrix
is therefore positive definite. Using (2), its determinant is

```
det [[m,p_1,p_2], [p_1,p_2,p_3], [p_2,p_3,p_4]]
  = (p_2^3/16)(1-rho)(9m rho-m-8) > 0.               (3)
```

Since `0<rho<1` (all `y_j` are nonzero and there is more than one),
(3) proves (1). Strictness follows from distinct roots, without an
asymptotic argument or an equality-case assumption.

## A separate bound from the paired sign order

The [critical root-order theorem](critical_odd_moment_root_order.md)
orders the signs, up to a common factor `epsilon`, as `++--++--...`.
Write `q_1>...>q_m>0` for their reciprocal absolute values. Grouping
terms in consecutive pairs makes the pair sums strictly decrease, so
the alternating sum obeys

```
0 < epsilon p_1 < q_1+q_2 < 2q_1.
```

For clarity, the final unpaired term causes no exception: the last
complete pair exceeds that final term, and all consecutive pair sums
decrease. Thus (2) also yields

```
p_2 < 4 max_j z_j^(-2).                              (4)
```

This improves the bound `p_2<9 max_j z_j^(-2)` that would follow from
(1) alone, but does not replace the fourth-moment assertion (1).
Neither inequality supplies a uniform reciprocal-tail decay estimate.

Using the entire odd-plus-constant residual identity, the later
[uniform-tail theorem](critical_moment_uniform_reciprocal_tail.md)
does prove degree-independent reciprocal-square tail tightness. Its
compactness proof gives no explicit rate. The subsequent
[harmonic proof](critical_moment_effective_reciprocal_tail.md) strengthens
this to linear magnitude growth and an explicit `O(1/K)` square tail.
Neither theorem transfers to arbitrary circle arcs.

The [exact checker](check_critical_moment_reciprocal_concentration.py)
verifies the determinant polynomial and a rational degree-five
critical fixture. The fixture is one support, not a full eight-row
solution.
