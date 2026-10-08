# Affine allocation rigidity and products of every length

Let distinct nonzero Gaussian integers `z_1,...,z_M` have modulus `R`
and lie on an arc of length at most `C sqrt(R)`. If `C<=sqrt(2)`,
their split-prime allocation vectors are **affinely independent over
the reals**, with arbitrary Gaussian units allowed. If the points have
one literal common Gaussian unit, the same assertion holds for `C<=2`.

Let `r` count only those rational split primes at which the allocation
varies among these points. Then

```text
M<=r+1.                                                (1)
```

For `R>1`, the points are multiplicatively independent even modulo
Gaussian units. In particular, on a `C<=sqrt(2)` arc, equality of two
products of any equal length implies equality of the multisets, with
repetitions allowed. The length has no upper bound and `C` need not
shrink as the length increases.

This strengthens the bounded-order conclusion of
[angular_moment_product_rigidity.md](angular_moment_product_rigidity.md)
and gives a linear conditional count in the number of varying split
primes. That number can grow with `R`, so (1) is not the requested
uniform point-count bound. Full balanced and fair patterns already
have full affine rank, as explained in
[relation_free_balanced_patterns.md](relation_free_balanced_patterns.md).

## 1. Exact pair separation with arbitrary units

Put `N=R^2` and `W=log N`. Fix one Gaussian prime `pi_p` above every
rational split prime dividing `N`, and write its allocation in row `i`
as `a_i(p)`, with `0<=a_i(p)<=v_p(N)`. The inert and ramified factors
are common. All remaining ambiguity is an explicit unit `epsilon_i`.

For two rows, the standard exact factorization gives

```text
z_i=epsilon_i g_ij A_ij,
z_j=epsilon_j g_ij bar(A_ij),
P_ij=Norm(A_ij)=product_p p^|a_i(p)-a_j(p)|,
|g_ij|^2=N/P_ij.                                      (2)
```

Write `A_ij=a+ib`, and put `eta=epsilon_j/epsilon_i`. Its chord
numerator `|A_ij-eta bar(A_ij)|` is one of

```text
2|b|,          2|a|,          sqrt(2)|a-b|,
sqrt(2)|a+b|,                                         (3)
```

for `eta=1,-1,i,-i`, respectively. The relevant integer is nonzero
because the two original points are distinct. Thus every pair obeys

```text
|z_i-z_j|^2 >= 2N/P_ij,
```

and in one common unit class the stronger lower bound is `4N/P_ij`.
The chord between distinct endpoints is strictly shorter than their
positive circular-arc gap, which is at most `C sqrt(R)`. Consequently

```text
log P_ij > W/2 + log(2/C^2)    for arbitrary units,
log P_ij > W/2 + log(4/C^2)    for one common unit.      (4)
```

At the stated constants, both give the strict inequality

```text
log P_ij>W/2.                                         (5)
```

The strictness at `C=sqrt(2)` or `C=2` is retained; no smaller constant
or large-radius threshold is being substituted.

## 2. Linear coordinates, with one coordinate per varying prime

For each varying split prime put

```text
m_p=min_i a_i(p),
e_p=max_i a_i(p)-m_p>0,
b_i(p)=a_i(p)-m_p,
W_0=sum_p e_p log p<=W.
```

The squared radius after removing the complete common Gaussian divisor
is `exp(W_0)`. All constant allocation coordinates have disappeared;
they are irrelevant to affine relations. Pair distances are unchanged:

```text
d_ij:=log P_ij=sum_p |b_i(p)-b_j(p)| log p.             (6)
```

For `0<=x,y<=e`, the elementary inequality

```text
|x-y| <= x+y-2xy/e                                    (7)
```

follows, when `x>=y`, by subtracting the left side to obtain
`2y(1-x/e)>=0`, and otherwise by symmetry.

One can therefore use the **linear allocation** vectors in `R^r`

```text
f_i(p)=sqrt(e_p log p) (2b_i(p)/e_p-1).
```

Unlike the earlier threshold-layer embedding, this uses only one
coordinate for each varying prime, regardless of its exponent. From
(7),

```text
<f_i,f_j> <= W_0-2d_ij <= W-2d_ij < 0    (i!=j).       (8)
```

Suppose the allocation rows had a nonzero real affine relation
`sum_i c_i b_i=0`, `sum_i c_i=0`. Because `b_i -> f_i` is an affine
map, `sum_i c_i f_i=0` as well. Split the nonzero coefficients into
positive and negative sets and write

```text
v=sum_(c_i>0)c_i f_i=sum_(c_j<0)(-c_j)f_j.
```

Both index sets are nonempty, and (8) gives

```text
||v||^2
 =sum_(c_i>0,c_j<0)c_i(-c_j)<f_i,f_j> < 0,
```

a contradiction. The allocation rows are affinely independent, proving
(1).

There is an equivalent convex formulation. Normalize the positive and
negative coefficients into disjoint probability measures with the same
allocation mean `mu_p`. Their mean cross-distance in coordinate `p`
is at most `2mu_p-2mu_p^2/e_p<=e_p/2` by (7). The mean of (6) is
therefore at most `W_0/2`, contradicting (5) for every cross pair.
This also proves the assertion for arbitrary real affine coefficients.

## 3. All-order product consequences and the unit-radius case

Assume `R>1` and suppose an integer relation has unit value:

```text
product_i z_i^n_i=u,             u in {1,-1,i,-i}.
```

Taking absolute values gives `sum_i n_i=0`. Taking the valuation at
each chosen split Gaussian prime gives `sum_i n_i a_i(p)=0`, hence
`sum_i n_i b_i(p)=0`. Affine independence forces every `n_i=0` and
then `u=1`. Thus the points are independent modulo Gaussian units,
which includes ordinary multiplicative independence.

The same argument immediately rules out unequal multisets in products
of any common length. This conclusion has no angular-sum unwrapping
step and no dependence on the product length.

When `R=1`, an arbitrary-unit arc of length at most `sqrt(2)` contains
at most one lattice point: distinct equal-norm points have chord at
least `sqrt(2)`, strictly shorter than their circular arc gap. A single
unit is torsion, so unrestricted multiplicative independence is false
in this case. Equal-length multiset rigidity is still immediate. A
single common-unit class likewise contains at most one point.

For any fixed larger `C`, subdividing into
`ceil(C/sqrt(2))` arcs gives the conditional count

```text
M <= ceil(C/sqrt(2)) (r+1),                             (9)
```

where `r` may be taken for the whole original tuple. Each subarc has
at most that many varying primes. Assign subdivision endpoints to one
piece to avoid double counting. Equation (9) does not assert all-order
product rigidity across different pieces.

## Verification scope

The proof is elementary and covers arbitrary prime powers, common
factors, real affine coefficients, and unit choices. The companion
[checker](check_linear_allocation_affine_rigidity.py) tests the exact
unit-dependent chord bounds and both allocation-rank conclusions
on finite actual Gaussian circles. It supplements the proof and is not
a proof by bounded search. No Lean formalization is claimed here.
