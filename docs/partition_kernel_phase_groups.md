# Several moving phase groups with a common norm budget

The rank-one argument in
[the relative-phase note](rank_one_kernel_relative_roth_obstruction.md)
extends to certain higher-dimensional kernels. The condition below is
imposed on the **full** row-difference space, so omitted blocks do not
silently become uncontrolled row phases. It also gives an explicit family
with arbitrarily many moving phase parameters and a uniform arithmetic
constant. No extraction from an arbitrary large endpoint cluster is proved.

## 1. A checkable condition on part of the full kernel

Use the fixed sign model

```text
Z_i = product_j kappa_j^((1+S_ij)/2)
                    bar(kappa_j)^((1-S_ij)/2),
z_i = d u_i Z_i,               R=|d| product_j |kappa_j|,
gcd_i Z_i=1 in Z[i].                                      (1)
```

Here d is a nonzero Gaussian integer, the u_i are arbitrary Gaussian
units, and the nonzero blocks kappa_j are coprime to their conjugates.
Primitivity in (1) concerns the whole normalized row set, not each row.
Set

```text
V=span_Q{S_i-S_0},                  K=V^perp.              (2)
```

Suppose disjoint coordinate sets I_g carry signs epsilon_j in {+1,-1}
such that, for every group,

```text
r_g=sum_(j in I_g) epsilon_j e_j lies in K,
epsilon_j e_j-epsilon_h e_h lies in V for j,h in I_g.       (3)
```

These are conditions in the entire k-coordinate space. In addition let
J be a set of coordinates with e_j in V. Such a coordinate cannot belong
to any I_g: its dot product with r_g would be nonzero.

Conjugate kappa_j and flip column S_j when epsilon_j=-1. Then (3) says
that each group has its own kernel indicator, and every coordinate
difference inside that group is a phase certificate. Each row has a
fixed plus-count q_g in I_g, since its scalar product with the group
indicator is fixed. The kernel outside these groups is unrestricted.

In particular (3) holds if, after sign changes, K is the space of vectors
constant on each set of a partition. It also holds for any partition
summand of a larger kernel; it is not necessary to control every column.

## 2. Primitive content supplies a matching in each group

First suppose every block of I_g is a nonunit. Partition these blocks
into associate classes after the preceding conjugations. A class of
size s satisfies

```text
s <= min(q_g, |I_g|-q_g) <= floor(|I_g|/2).              (4)
```

Choose a Gaussian prime pi dividing a class member. It is not associated
to its conjugate, and no class member contains bar(pi). Whole-set
primitivity supplies one pi-free row and one bar(pi)-free row. They must
choose, respectively, all minus signs and all plus signs on the class.
Other blocks have nonnegative valuations and cannot cancel this
requirement. The fixed plus-count proves (4), including when other
groups share the same rational prime.

Place associate classes contiguously on a cyclic list of n=|I_g| vertices.
For odd n use edges i--i+(n-1)/2 modulo n, each with weight 1/2. For even
n use antipodal pairs, each with weight one. Condition (4) makes every
edge a nonassociate pair. Each vertex has total edge weight one, and
the total edge weight is n/2.

For each edge jh, (3) removes the moving phase parameter. The directional
reduction of kappa_j bar(kappa_h) is a nonunit exactly because the two
blocks are nonassociate. Its logarithmic norm is at most
log Norm(kappa_j)+log Norm(kappa_h). Exact shared-prime cancellation only
reduces this upper bound. The fixed rational certificate gives a finite
set of rational-angle targets with error O_S(Delta), where Delta is the
angular span of the original points.

If I_g contains a unit coordinate h, all its nonunit coordinates instead
become individually isolable: use the pair certificate with h and absorb
the unit angle into the finite target set. Each such coordinate has cost
one. Formally one can delete all unit columns and adjust the row units;
the resulting individual certificate is then a standard coordinate vector.
This avoids charging arbitrary capacity to a zero-height unit block.

## 3. General group theorem and its growth exponent

Let A count the nonunit coordinates in J and in groups that contain at
least one unit. Let B count the coordinates in groups whose blocks are
all nonunits. All these sets are disjoint. Define

```text
P=A+B/2.                                                (5)
```

For each fixed S satisfying (3), fixed C>0, and fixed pattern of unit
coordinates, the conditions (1) and P>4 exclude unbounded R on arcs
of length at most C sqrt(R). In particular nine nonunit coordinates
covered by these groups or by J suffice. Unit patterns and associate
partitions have only finitely many possibilities for fixed S, so the
conclusion remains uniform over them.

More precisely, for every eta>0 there is c=c(S,eta)>0 such that

```text
L=R Delta >= c R^(1-(2+eta)/P).                         (6)
```

To prove this, fixed-target Roth separation gives, for each counted
nonunit reduced product g,

```text
log Norm(g) >= [2/(2+eta)] log(1/Delta)-O_(S,eta)(1).
```

Sum with the matching and isolation weights. Every nonunit block has
load at most one, so the sum is at most 2 log R. The total weight is P,
which proves (6). Large Delta can be included by reducing c. All
constants in this general argument may be ineffective. This invokes
only the fixed-target Roth lemma already stated and sourced in
[the coupled-product theorem](coupled_block_roth_packing.md); it is a
consequence of that packing theorem with a proved nonunit criterion.

For P>4 choose eta<P/2-2. The exponent in (6) then exceeds 1/2.
The groups can have several independent moving parameters: no one of
those parameters is treated as a fixed algebraic target.

## 4. An effective family with arbitrarily many moving parameters

Let T be the 11-by-11 core of a normalized order-twelve Hadamard matrix.
Use b>=1 groups of eleven nonunit conjugate-primitive blocks. The
following sign matrix has 10b+1 rows and 11b columns.

Start with the state (a_1,...,a_b)=(0,...,0). For h=1,...,10, increment
each of a_1,...,a_b in that order from h-1 to h, recording the state after
every increment. At each recorded state concatenate core rows

```text
S_state=(T_(a_1),...,T_(a_b)).                           (7)
```

Successive row differences change just one group from T_(h-1) to T_h.
Those ten differences span the ten-dimensional sum-zero space in that
group. Consequently

```text
dim K=b,
K=span_Q{the b group indicators}.                       (8)
```

For b=2 this is a 21-row, 22-column profile with two moving phases.
The interleaving matters to the example: when b>1, the rows are not a
Cartesian product of full cores with all other groups frozen.

For coordinates a!=c in a group, put

```text
l_i=(T_(i,a)-T_(i,c))/12,       i=0,...,10.
```

The core identities give sum l_i=0, sum |l_i|=1, and
sum_i l_i T_i=e_a-e_c. Let y_t denote a lifted full-row phase. For each
h let t(g,h) be the row index immediately after incrementing group g
from h-1 to h. Its immediately preceding row has index t(g,h)-1.
The full-space certificate is

```text
sum_(h=1)^10 (sum_(i=h)^10 l_i)
              (y_(t(g,h))-y_(t(g,h)-1)).                (9)
```

Indeed these differences remove every other group. Telescoping within
group g gives sum_i l_i T_i theta_g. Every coefficient in (9) is an
integer multiple of 1/6, the total coefficient sum is zero, and the
l1 norm after collecting coefficients is at most

```text
2 sum_(h=1)^10 |sum_(i=h)^10 l_i|
 <= 2 sum_(i=1)^10 i |l_i| <= 20.                      (10)
```

Neither the denominator nor this bound depends on b. Thus integer phase
lifts and arbitrary independent row units give the exact target set
(pi/12)Z, and the pair phase error is at most 10 Delta.

The elementary fifteen-degree separation in the relative-phase note is

```text
dist(arg g,(pi/12)Z) >= 1/(4 Norm(g))
```

for every conjugate-primitive nonunit g. All core rows have five plus
signs, so (4) bounds associate classes in each group by five. Use the
eleven-cycle matching separately in every group. Every counted pair
therefore satisfies

```text
Norm(kappa_j) Norm(kappa_h) >= 1/(40 Delta).             (11)
```

There are b groups, each with total matching weight 11/2. The common
norm budget from (1), which retains a common multiplier d, now gives

```text
R^2 >= product_j Norm(kappa_j)
    >= (1/(40 Delta))^(11b/2).
```

Hence the explicit family (7), under whole-set primitivity, satisfies

```text
L >= (1/40) R^(1-4/(11b)).                              (12)
```

For b=2 the growth exponent is 9/11. For endpoint arcs it follows that

```text
R <= (40C)^(22b/(11b-8)).                              (13)
```

Every conjugate-primitive nonunit has norm at least five. Combining this
with (12) also bounds the number of groups directly:

```text
C >= (1/40) 5^((11b-8)/4).                             (14)
```

Thus the constants here are uniform even as this particular family of
matrices varies with b. Shared Gaussian primes, repeated block directions,
independent row units, and a common Gaussian multiplier are all permitted
under (1). Equations (12)--(14) do not assert that every large cluster
contains such a profile or can be regrouped into one.

The accompanying exact checker verifies the full-space certificates,
their denominators and loads, and the higher-dimensional kernels for
several b; the proof above covers every b.
