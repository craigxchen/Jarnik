# Auxiliary forms and the dual lattice of the Gaussian metric lift

This note investigates the rank-one Gaussian metric lift beyond its
individual Pluecker equations. It establishes an exact description of
the private-factor divisibility of short relation minors, and a dimension
and volume barrier for using those minors to prove the endpoint bound.
No uniform cardinality theorem is proved.

The starting metric reconstruction is the one in Section 6.8 of
`codex_uniform_bound_research.md`: Gaussian integers `B_1,...,B_m`
lie on nearby rays, their norms have specified private factors, and
the Hermitian matrix is exactly `(conjugate(B_i) B_j)`. All auxiliary
linear relations below use both coordinates of these genuine Gaussian
integers; they are not relations in an abstract Pluecker relaxation.

## 1. A first small-value calculation lands exactly at the fair-model scale

In the complete fair model, `m>=5` is odd, each private Gaussian block
occurs in exactly `t=(m-1)/2` of the `B_i`, and the normalized
logarithmic block weights tend to equality. The exact fractions below
refer to this limiting exponent model; equal norms of distinct coprime
Gaussian blocks are not assumed. Then

```text
|B_i| = R^(alpha+o(1)),       alpha=(m-1)/(2m),
angular width = O_C(R^(-1/2)).
```

Choose axes along their common ray. The longitudinal coordinates have
size `R^(alpha+o(1))` and the transverse coordinates have size
`R^(alpha-1/2+o(1))`. Pigeonholing bounded integer combinations of
the longitudinal coordinates gives a nonzero coefficient vector of
size `O(K)` with longitudinal sum

```text
R^(alpha+o(1))/K^(m-1).
```

Its transverse sum has the independent bound

```text
K R^(alpha-1/2+o(1)).
```

To make both decay as a power of `R` would require simultaneously

```text
K > R^(alpha/(m-1)) = R^(1/(2m)),
K < R^(1/2-alpha)   = R^(1/(2m)).
```

There is no strict power interval. A small enough arc constant can
still produce exact Gaussian linear relations at the boundary scale.
That possibility must be retained: it is insufficient simply to say
that `m>2` Gaussian integers already admit rational linear relations.
The next sections check the sizes and the private-factor constraints
of several independent relations.

## 2. The exact relation lattice and its covolume

Let `V` be the integer `2 by m` matrix whose columns are
`(Re B_i,Im B_i)`. The nearby rays are distinct, so `V` has rank two.
Set

```text
Lambda = {c in Z^m : sum_i c_i B_i=0},
r = m-2,
D_ij = det(B_i,B_j),
g = gcd_(i<j) |D_ij|.
```

The relation lattice has rank `r`, is saturated in its rational span,
and has Euclidean covolume

```text
det(Lambda) = sqrt(sum_(i<j) D_ij^2)/g.             (1)
```

Equivalently, if `K` is the `r by m` matrix of an integer basis of
`Lambda`, its complementary maximal minors are

```text
det K_([m] minus {i,j}) = +/- D_ij/g.              (2)
```

One way to see this is to take exterior coordinates of the rational
kernel: orthogonal complementation exchanges the two-row minors with
the complementary kernel minors. Saturation makes the latter primitive,
so the common divisor `g` is exactly the scalar removed. Cauchy--Binet
then proves (1).

The product of the successive minima is bounded above and below by
dimension-dependent multiples of (1). Thus the metric lift controls a
whole product of short relation sizes, not just one relation.

## 3. Private supports force lower-order coefficient minors as well

Suppose a private block has norm `p^e`, where `p` is a split rational
prime, and occurs in exactly the nonempty proper set `T` of columns.
Choose its Gaussian orientation as `pi^e`. The block assumptions give

```text
pi^e divides B_i exactly when i belongs to T,
conjugate(pi) does not divide any B_i,
pi does not divide B_i when i is outside T.
```

In particular the two-coordinate matrix has rank two modulo `p`:
an inside column is a nonzero vector on one isotropic line, while an
outside column is not on that line. Therefore `p` does not divide `g`
in (1). In the complete fair model all conductor primes are active in
this sense, so `g` is coprime to the conductor.

Now take any `h` relation vectors in `Lambda`, and any `h`-element
column set `J`. If

```text
J^c is contained in T,                             (3)
```

then `p^e` divides the corresponding `h by h` coefficient minor.
To prove this, reduce the relation `sum c_i B_i=0` modulo `pi^e`.
Since `Z[i]/(pi^e)` is `Z/p^e Z`, this is one linear congruence on
the coefficients. Its entries vanish on `T`, and are units outside
`T`. Condition (3) puts its entire nonzero support inside `J`, giving
a nonzero, unimodular vector in the right kernel of the square
coefficient submatrix modulo `p^e`. Elementary column operations over
`Z/p^e Z` then show that its determinant is zero modulo `p^e`.

Consequently every such coefficient minor is divisible by

```text
Q_J = product_(private supports T with J^c subset T) N(H_T).   (4)
```

This applies to any independent short relations, whether or not they
are already a basis. Since some `h`-minor is nonzero, Hadamard gives
the useful inequality

```text
product_(l=1..h) ||c_l|| >= min_(|J|=h) Q_J.        (5)
```

The order `h=m-2` is exterior duality applied to the old pair-minor
divisibility. Lower orders are additional constraints on an individual
short basis, so they must be checked before declaring the route closed.

## 4. In the complete fair model all these constraints remain feasible

Let the normalized logarithmic weights tend to the uniform distribution
on the `t`-subsets of `[m]`, where `t=(m-1)/2`. The total log norm of
all blocks is `2 log R`. For a
fixed set of `u` rows, the fraction of supports containing it is

```text
a_u = binom(m-u,t-u)/binom(m,t) = (t)_u/(m)_u,
```

with `a_u=0` when `u>t`. Thus every divisor in (4), for `|J|=h`,
has the limiting size

```text
Q_J = R^(2 a_(m-h)+o(1)).                          (6)
```

The metric upper bound for a pair minor is

```text
|D_ij| <= R^(beta+o(1)),
beta=2alpha-1/2=(m-2)/(2m).
```

Equations (1) and the successive-minima product therefore permit the
balanced scale

```text
||c_l|| = R^(1/(2m)+o(1))       for every l=1,...,m-2,
```

when `g=R^o(1)`. This is a feasible exponent assignment, not an
assertion that every actual lattice has equal successive minima.

At that assignment, every lower-minor constraint (5) is consistent,
because

```text
2 a_(m-h) < h/(2m)       for 1<=h<=m-2.            (7)
```

For a direct proof, put `u=m-h>=2`. At `u=2`,

```text
a_2=(m-3)/(4m) < (m-2)/(4m).
```

The recurrence `a_(u+1)/a_u=(t-u)/(m-u)` preserves
`a_u<(m-u)/(4m)` until the numerator becomes zero. This proves (7).

The strongest constraint is the maximal-minor one. Its divisor and
permitted size are respectively

```text
R^((m-3)/(2m))       and       R^((m-2)/(2m)).
```

The gap is exactly `R^(1/(2m))`, the existing primitive-residue scale.
In particular, using many short relations and then taking their
maximal minors reproduces the old cofactor gap through (2).

Any positive common divisor `g` reduces the permitted product of
successive minima. The maximal-minor constraint then gives at most
`g <= R^(1/(2m)+o(1))`. It does not force `g` to grow, and the case
`g=1` remains compatible with every displayed exponent constraint.

Thus no nonnegative combination or optimization of these local minor
lower bounds and the covolume upper bound yields a strict exponent
gain: the common feasible exponent assignment above witnesses the
barrier. This does not assert existence of a Gaussian endpoint family
with all these data, nor exclude a different use of residue values.

## 5. The uniform-profile model has zero leading margin

The following calculation retains the improvement from the tuple's
own common factor.

There is a sharper comparison with the current extraction target.
Here the metric block calculation is made in the squarefree or
endpoint-allocation model, where private supports are independent.
It must not be applied automatically to arbitrary prime-power
threshold layers: distinct layers can involve the same prime.

Select `k=m+1` rows with inherited asymptotically uniform pattern
weights and divide their own common Gaussian gcd. Orient every cut
away from the anchor row. For `B_i=G/A_i`, the private supports `T`
now range over all **proper** subsets of `[m]`, including the empty
set, each with intrinsic weight fraction `1/(2^m-1)`. Write the
intrinsic radius as `R'` and put `d=2^m-1`.

The removed inherited weight fraction is `2^(-m)+o(1)`. Keeping the
resulting improvement in the arc constant gives

```text
|B_i| = R'^(alpha+o(1)),       alpha=(2^(m-1)-1)/d,
Delta <= R'^(-rho+o(1)),       rho=2^(m-1)/d.
```

Thus the pair-minor exponent is

```text
beta=2alpha-rho=2(2^(m-2)-1)/d.                   (8)
```

For an `h`-element set `J`, the fraction of proper supports containing
`J^c` is `(2^h-1)/d`. Accordingly the coefficient-minor divisor in
(4) has exponent

```text
2(2^h-1)/d.                                      (9)
```

The balanced successive-minimum assignment is now

```text
log ||c_l||/log R' = beta/(m-2)+o(1).
```

It satisfies every local minor constraint, since

```text
(2^h-1)/h <= (2^(m-2)-1)/(m-2),
```

with equality exactly when `h=m-2`. The ratio is strictly increasing
for positive integer `h`: the difference of successive numerators
after cross-multiplication is `(h-1)2^h+1>0`.

Thus the maximal-minor condition **saturates** the leading covolume
exponent; the smaller-minor conditions have slack. This is zero
leading margin, rather than the fixed positive gap in Section 4.
The `R'^o(1)` cofactors allowed by extraction can still absorb fixed
constants, so leading equality is not a contradiction.

The empty support divides none of the `B_i`; its primes need not be
coprime to `g`. All nonempty proper supports are active, and the
argument from Section 3 still applies to them. The maximal-minor
comparison forces only `g=R'^o(1)` in this limiting model, which is
again compatible with `g=1`.

## 6. Why simple positivity or a modular quadratic does not fix this

An auxiliary integer form can be chosen in a congruence lattice so
that its value is divisible by a large conductor. A volume argument
can then find signed coefficients making its archimedean value small.
For the linear metric lift, the vanishing obtained in this manner
belongs to the relation lattice analyzed above. Its coefficients and
their minors cannot be treated as bounded independently of the radius.

There is a separate obstruction to adding positivity without proof.
For example, let `q_ij=|z_i-z_j|^2>0` and take nonnegative integer
coefficients `a_ij`, not all zero. If

```text
F=sum_(i<j) a_ij q_ij
```

is divisible by the squared radius `N=R^2`, then automatically
`F>=N`. On an endpoint arc, `q_ij<=C^2 R`, so any such congruence
solution satisfies

```text
sum_(i<j) a_ij >= R/C^2.
```

Minkowski's theorem for a symmetric coefficient body cannot supply
a positive coefficient vector smaller than this. Replacing its
signed small vector by a positive one would be an additional theorem,
and here the elementary value bound shows precisely why it fails.

The current output is therefore a checked limitation of these specified
auxiliary-form operations. A successful new argument would need to use
joint residue equations, a proved restriction on successive-minimum
shapes beyond the displayed minors, or an auxiliary expression whose
vanishing has a stronger consequence than these existing relations.

## Verification

An exact seven-row Gaussian family with one distinct split prime on
each three-element private support verified 196 coefficient-minor
divisibilities from (4), using several independent integer relations
and both nontrivial minor orders. Exact rational checks verified (7)
and the uniform-profile ratio inequality for every odd `m` from 5
through 101. These supplement the general proofs; the test family
is not asserted to lie on an endpoint arc.
