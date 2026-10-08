# One spectral coefficient orders every critical pair orientation

The critical and successor root-order rules have a further joint
constraint: their first-sign orientations must come from one strict
ordering of eight real spectral coefficients. This is valid in every
degree in the punctured-orthogonal setup. It does not by itself exclude
that setup; a simultaneous ordering or path obstruction is still needed.

Let `T` have eight orthogonal sign rows and `4s+4` columns. Delete a
constant column and a balanced column, leaving `S` with `n=4s+2`
columns. Let `A,B` denote the two four-row groups defined by the deleted
balanced column. Suppose nonzero real `a_j` with pairwise distinct
absolute values obey

```
S(a^(2r-1)) = alpha_r 1,       1<=r<=s.
```

Write `d=2s+1` and define

```
P_i(t) = product_(j=1)^n (t-S_ij a_j),
kappa_i = coefficient of t^d in P_i(t).
```

All the polynomials `P_i` have the same first `2s=d-1` coefficients
below their monic leading term. Their even power sums coincide because
the squared roots are independent of the row; their odd power sums
through degree `2s-1` coincide by hypothesis. Newton's identities give
the coefficient assertion. In addition,

```
P_i(t) P_i(-t) = product_j (t^2-a_j^2)                 (1)
```

for every row, since `n` is even.

## Literal differences of spectral factors

For rows `i,j`, let `D_ij` be their differing-column set and set

```
G_ij(t)=product_(k notin D_ij)(t-S_ik a_k),
R_ij(t)=product_(k in D_ij)(t-S_ik a_k).
```

Cross-group pairs have `|D_ij|=d`; within-group pairs have
`|D_ij|=d+1`. The roots of `R_ij` have zero odd power sums through
degree `2s-1`. The critical and successor Newton reductions therefore
give respectively

```
cross:  R_ij(t)=F_ij(t)+r_ij,         F_ij odd,
within: R_ij(t)=Q_ij(t^2)+r_ij t.
```

Here `r_ij!=0`: in the critical case it is the nonzero constant
coefficient; in the successor case zero would give opposite roots,
contrary to distinct absolute values. Consequently

```
cross:  P_i-P_j = 2r_ij G_ij,
within: P_i-P_j = 2r_ij t G_ij,
kappa_i-kappa_j = 2r_ij != 0.                         (2)
```

Thus all eight `kappa_i` are distinct. Every pair difference is a
real-rooted polynomial of exact degree `d`, with its roots precisely
the common signed source roots, and also zero for a within-group pair.
No content normalization or asymptotic approximation is involved.

## The first signs determine a transitive ordering

Orient a certificate by `c=(S_i-S_j)/2`. Let `epsilon_ij` be the
first sign of `c_k a_k` when its support is read by increasing
absolute value. The proved
[critical and successor lemma](critical_odd_moment_root_order.md)
then gives

```
pair type       sign(kappa_i-kappa_j)
critical        (-1)^(s+1) epsilon_ij
successor A     (-1)^s     epsilon_ij
successor B     (-1)^(s+1) epsilon_ij.                  (3)
```

For the critical case the central root in the odd-fiber proof has
sign `-sign(r_ij)(-1)^s`. For successor A, `Q_ij` has `s+1` positive
roots, so `sign Q_ij(0)=(-1)^(s+1)`; the smallest root of
`Q_ij(t^2)+r_ij t` has sign `-sign(Q_ij(0)/r_ij)`.
For successor B, `Q_ij` instead has `s` positive roots and one negative
root, reversing that constant-coefficient sign. These facts give (3),
and (2) transfers the sign of `r_ij` to the coefficient difference.

In particular the tournament obtained from the right-hand signs in
(3) must be transitive. A directed cycle is a contradiction, regardless
of the degree, coefficient magnitudes, or prime allocation model.
Independent admissible pair patterns need not satisfy this condition.

## Reduction after fixing the final cross orientations

Put `u_j=sign(a_j)` and `H_i=sum_j S_ij u_j`. A critical pattern has
total sign `epsilon_ij`, so for cross pairs

```
(H_a-H_b)/2=epsilon_ab in {+1,-1}.
```

Every cross-group distance between these scalar heights is therefore
two. At least one group has a common height: if two heights in `A`
differ, their common distance-two neighbors in `B` consist of their
unique midpoint. After interchanging groups, assume all `H_a`, `a in A`,
are equal. Then `epsilon_ab=epsilon_b` depends only on `b`.

Within `A`, and within either fixed-`epsilon_b` subgroup of `B`, the
total certificate sign is zero, so the pair must have successor type A.
Between opposite-`epsilon_b` subgroups the total is two in absolute
value, so the pair must have successor type B. Thus all pair types are
determined by the final cross orientations.

Equation (3) places the four coefficients `kappa_a` in one common
interval between the two `epsilon_b` classes. More explicitly,

```
sign(kappa_a-kappa_b)=(-1)^(s+1) epsilon_b
```

for every `a`. Internal orders in these three groups are unrestricted
by this statement. If the two classes of `B` have sizes `q` and `4-q`,
there are at most

```
4! q! (4-q)!
```

possible full coefficient orders for fixed labeled cross orientations.
Each such order fixes all successor-A first signs by (3). This replaces
independent pair-orientation choices by a bounded family of transitive
orders, uniformly in `s`.

## Initial directed cuts

Direct a pair `i -> j` when its required first certificate sign
`epsilon_ij` is positive. This initial tournament is itself transitive:
between the three coefficient blocks its direction agrees with
`(-1)^(s+1)(kappa_i-kappa_j)`, while inside each block it agrees with
`(-1)^s(kappa_i-kappa_j)`. Thus order the blocks one way and their
internal coefficients the other way. Empty blocks cause no exception.

For the first column after absorbing its coefficient sign, let
`v_i=S_ij sign(a_j)`. Whenever `v_i=+1,v_k=-1`, its first certificate
sign is positive, so every edge across this cut must point from the
plus set to the minus set. The plus set must therefore be an initial
segment of the initial tournament order. There are only seven possible
nonconstant first signed columns for a fixed order.

More generally, at a later prefix orient each pair by the next required
sign in its critical or successor pattern. A permissible signed column
is exactly a directed cut of this tournament. Its plus set is an initial
union of strongly connected components, whose condensation is totally
ordered. In particular a strongly connected next-sign tournament permits
no nonconstant next column. This is an exact finite-state reduction,
not an assumption of global transitivity after the first column.

The corresponding prefix-height constraints remain necessary: the four
`A` heights, divided by two, have diameter at most one; each fixed-sign
`B` class also has diameter at most one; and the opposite `B` classes
are separated by a nonnegative height difference of at most three.
They follow from the partial sums of successor A and B. A finite-state
exclusion retaining these constraints and the individual pair distances
would be additional work. Neither the coefficient order nor these
bounded prefix heights alone establish the desired general exclusion.

[The exact checker](check_critical_moment_spectral_orientation.py)
checks residual-polynomial differences and coefficient signs in rational
real-rooted critical and successor fixtures. It also checks the initial
tournament and all 256 candidate signed columns for both parities of
`s` and every division of the four `B` rows. These are pair and ordering
fixtures; they do not assert existence of a complete eight-row moment
solution. The proof above requires no Hadamard completion or additional
genericity of the coefficients `kappa_i`.

## A transitive value ordering at every prefix

There is a stronger prefix test separate from the next-sign tournament.
Put `v_ij=S_ij sign(a_j)` and order the columns by increasing `|a_j|`.
At a prefix, let `u_i` count processed positive entries in row `i`,
and let `k_ij` count processed differing entries for the pair. Let
`p_i=(n+H_i)/2` be the total positive-entry count, and let `D_ij` be
the full pair distance. The number of unprocessed common positive roots is

```
M_ij = (p_i+p_j-D_ij-u_i-u_j+k_ij)/2.                  (4)
```

For a positive `t` between this prefix's last magnitude and the next
magnitude, the exact factorization (2) gives

```
sign(P_i(t)-P_j(t))
 = sign(kappa_i-kappa_j) (-1)^M_ij.                    (5)
```

Only unprocessed common positive roots contribute negative factors in
`G_ij(t)`; the extra factor `t` for within-group pairs is positive.
Thus the tournament prescribed by (5) must be transitive at every
prefix. None of the values `P_i(t)` coincide there, since all pair
differences have only the displayed shared roots and possibly zero.

The unknown common final height does not obstruct this test. With the
normalization above, write `p_i=p_0+delta_i`, where `delta_a=0` on `A`
and `delta_b=-epsilon_b` on `B`. Changing `p_0` reverses all signs in
(5) together when its parity changes, which preserves transitivity.
The integers in (4), up to this common shift, are therefore determined
by the prefix and the fixed cross orientations. This is a necessary
test on every prefix of an actual moment solution; it is not a claim
that a path passing the tests is realizable by polynomial roots.
